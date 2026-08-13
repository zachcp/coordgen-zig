const std = @import("std");
const core = @import("core");
const model = @import("model");
const optimize = @import("optimize");

const crown_grid_interval: f32 = 20;
const first_crown_distance: f32 = 60;
const crown_spacing: f32 = 60;
const pocket_border_allowance: f32 = 10;

/// Build upstream's `crown_index`-th ordered contour around the ligand. The
/// returned largest contour is allocator-owned. Pocket-distance metadata is
/// not part of the stable Zig input yet, so callers supply one value per atom;
/// zeroes reproduce ordinary ligand behavior.
pub fn shapeAroundLigand(
    allocator: std.mem.Allocator,
    atoms: []const core.math.Vec2,
    bonds: []const LigandBond,
    pocket_distances: []const f32,
    crown_index: u32,
) core.errors.Error![]core.math.Vec2 {
    if (atoms.len == 0 or atoms.len != pocket_distances.len) return error.InvalidMapping;

    var min = atoms[0];
    var max = atoms[0];
    var max_pocket_distance: f32 = 0;
    for (atoms, pocket_distances) |atom, pocket_distance| {
        if (!atom.isFinite() or !std.math.isFinite(pocket_distance) or pocket_distance < 0) {
            return error.InvalidCoordinate;
        }
        min.x = @min(min.x, atom.x);
        min.y = @min(min.y, atom.y);
        max.x = @max(max.x, atom.x);
        max.y = @max(max.y, atom.y);
        max_pocket_distance = @max(max_pocket_distance, pocket_distance);
    }
    for (bonds) |bond| {
        if (bond.start >= atoms.len or bond.end >= atoms.len) return error.InvalidMapping;
    }

    const border = crown_spacing * @as(f32, @floatFromInt(crown_index)) + first_crown_distance;
    if (!std.math.isFinite(border)) return error.TooManyItems;
    const padding = border + max_pocket_distance + pocket_border_allowance;
    const origin = core.math.Vec2{ .x = min.x - padding, .y = min.y - padding };
    const upper = core.math.Vec2{ .x = max.x + padding, .y = max.y + padding };
    if (!origin.isFinite() or !upper.isFinite()) return error.InvalidCoordinate;

    const width = try gridDimension(upper.x - origin.x);
    const height = try gridDimension(upper.y - origin.y);
    const value_count = std.math.mul(usize, width, height) catch return error.TooManyItems;
    const values = allocator.alloc(f32, value_count) catch return error.OutOfMemory;
    defer allocator.free(values);

    for (0..height) |y| {
        for (0..width) |x| {
            const point = core.math.Vec2{
                .x = origin.x + @as(f32, @floatFromInt(x)) * crown_grid_interval,
                .y = origin.y + @as(f32, @floatFromInt(y)) * crown_grid_interval,
            };
            var shortest = std.math.inf(f32);
            for (atoms, pocket_distances) |atom, pocket_distance| {
                shortest = @min(shortest, distance(atom, point) - pocket_distance - border);
            }
            for (bonds) |bond| {
                const segment = pointSegmentDistance(point, atoms[bond.start], atoms[bond.end]);
                const interpolated_pocket_distance = pocket_distances[bond.start] * (1 - segment.parameter) +
                    pocket_distances[bond.end] * segment.parameter;
                shortest = @min(shortest, segment.distance - interpolated_pocket_distance - border);
            }
            values[x + y * width] = shortest;
        }
    }

    var contour = try optimize.marching_squares.run(
        allocator,
        values,
        width,
        height,
        origin,
        .{ .x = crown_grid_interval, .y = crown_grid_interval },
        0,
    );
    defer contour.deinit();
    var ordered = try contour.orderedCoordinates(allocator);
    defer ordered.deinit();
    if (ordered.count() == 0) return allocator.alloc(core.math.Vec2, 0) catch return error.OutOfMemory;

    var largest_index: usize = 0;
    for (1..ordered.count()) |index| {
        if (ordered.contour(index).len > ordered.contour(largest_index).len) largest_index = index;
    }
    return allocator.dupe(core.math.Vec2, ordered.contour(largest_index)) catch return error.OutOfMemory;
}

pub const LigandBond = struct {
    start: u32,
    end: u32,
};

/// Place all residue representative atoms in successive crowns around an
/// already-laid-out ligand. Protein-only inputs use a separate chain-grid
/// phase and are rejected here.
pub fn placeAroundLigand(
    allocator: std.mem.Allocator,
    atoms: []model.Atom,
    bonds: []const model.Bond,
    residues: []const model.Residue,
    chain_bytes: []const u8,
    interactions: []const model.ResidueInteraction,
) core.errors.Error!void {
    if (residues.len == 0) return;
    const residue_atoms = allocator.alloc(bool, atoms.len) catch return error.OutOfMemory;
    defer allocator.free(residue_atoms);
    @memset(residue_atoms, false);
    for (residues, 0..) |residue, index| {
        if (residue.id.index() != index or residue.atom.index() >= atoms.len or residue_atoms[residue.atom.index()]) {
            return error.InvalidMapping;
        }
        residue_atoms[residue.atom.index()] = true;
    }

    var ligand_ids: std.ArrayList(core.ids.AtomId) = .empty;
    defer ligand_ids.deinit(allocator);
    var ligand_coordinates: std.ArrayList(core.math.Vec2) = .empty;
    defer ligand_coordinates.deinit(allocator);
    const ligand_index = allocator.alloc(u32, atoms.len) catch return error.OutOfMemory;
    defer allocator.free(ligand_index);
    @memset(ligand_index, std.math.maxInt(u32));
    for (atoms, residue_atoms, 0..) |atom, is_residue, index| {
        if (is_residue) continue;
        if (ligand_ids.items.len >= std.math.maxInt(u32)) return error.TooManyItems;
        ligand_index[index] = @intCast(ligand_ids.items.len);
        ligand_ids.append(allocator, core.ids.AtomId.fromIndex(@intCast(index))) catch return error.OutOfMemory;
        ligand_coordinates.append(allocator, atom.coordinates) catch return error.OutOfMemory;
    }
    if (ligand_ids.items.len == 0) return error.Unsupported;

    var ligand_bonds: std.ArrayList(LigandBond) = .empty;
    defer ligand_bonds.deinit(allocator);
    for (bonds) |bond| {
        if (bond.start.index() >= atoms.len or bond.end.index() >= atoms.len) return error.InvalidMapping;
        const start = ligand_index[bond.start.index()];
        const end = ligand_index[bond.end.index()];
        if (start != std.math.maxInt(u32) and end != std.math.maxInt(u32)) {
            ligand_bonds.append(allocator, .{ .start = start, .end = end }) catch return error.OutOfMemory;
        }
    }
    const pocket_distances = allocator.alloc(f32, ligand_ids.items.len) catch return error.OutOfMemory;
    defer allocator.free(pocket_distances);
    @memset(pocket_distances, 0);

    var groups = try groupSecondaryStructures(allocator, residues, chain_bytes);
    defer groups.deinit();
    const group_order = allocator.alloc(u32, groups.count()) catch return error.OutOfMemory;
    defer allocator.free(group_order);
    for (group_order, 0..) |*group, index| group.* = @intCast(index);
    std.mem.sortUnstable(u32, group_order, GroupPriorityContext{
        .groups = groups,
        .residues = residues,
        .interactions = interactions,
    }, groupPriorityLessThan);

    const coordinates_set = allocator.alloc(bool, atoms.len) catch return error.OutOfMemory;
    defer allocator.free(coordinates_set);
    for (coordinates_set, residue_atoms) |*is_set, is_residue| is_set.* = !is_residue;

    var remaining = residues.len;
    var crown_index: u32 = 0;
    while (remaining != 0) {
        const shape = try shapeAroundLigand(allocator, ligand_coordinates.items, ligand_bonds.items, pocket_distances, crown_index);
        defer allocator.free(shape);
        if (shape.len == 0) return error.InvalidMapping;
        const penalties = allocator.alloc(bool, shape.len) catch return error.OutOfMemory;
        defer allocator.free(penalties);
        @memset(penalties, false);
        const before = remaining;
        remaining = 0;
        for (group_order) |group_index| {
            remaining += try placeSecondaryStructure(groups.members(group_index), .{
                .allocator = allocator,
                .atoms = atoms,
                .residues = residues,
                .interactions = interactions,
                .ligand_atoms = ligand_ids.items,
                .coordinates_set = coordinates_set,
                .shape = shape,
                .shape_number = crown_index + 1,
                .penalties = penalties,
            });
        }
        if (remaining >= before) return error.InvalidMapping;
        crown_index = std.math.add(u32, crown_index, 1) catch return error.TooManyItems;
    }
}

const GroupPriorityContext = struct {
    groups: SecondaryStructures,
    residues: []const model.Residue,
    interactions: []const model.ResidueInteraction,
};

fn groupPriorityLessThan(context: GroupPriorityContext, first: u32, second: u32) bool {
    const first_score = groupPriority(context, first);
    const second_score = groupPriority(context, second);
    return if (first_score != second_score) first_score > second_score else first < second;
}

fn groupPriority(context: GroupPriorityContext, group_index: u32) f32 {
    const members = context.groups.members(group_index);
    var interaction_count: usize = 0;
    for (members) |id| {
        const atom = context.residues[id.index()].atom;
        for (context.interactions) |interaction| {
            interaction_count += @intFromBool(interaction.start == atom or interaction.end == atom);
        }
    }
    const length: f32 = @floatFromInt(members.len);
    return length + 3 * @as(f32, @floatFromInt(interaction_count)) / length;
}

fn gridDimension(extent: f32) core.errors.Error!usize {
    if (!std.math.isFinite(extent) or extent <= 0) return error.InvalidCoordinate;
    const cells = @floor(extent / crown_grid_interval);
    if (cells > @as(f32, @floatFromInt(std.math.maxInt(usize) - 2))) return error.TooManyItems;
    return @as(usize, @intFromFloat(cells)) + 2;
}

fn distance(first: core.math.Vec2, second: core.math.Vec2) f32 {
    return @sqrt(distanceSquared(first, second));
}

fn distanceSquared(first: core.math.Vec2, second: core.math.Vec2) f32 {
    const dx = second.x - first.x;
    const dy = second.y - first.y;
    return dx * dx + dy * dy;
}

const SegmentDistance = struct { distance: f32, parameter: f32 };

fn pointSegmentDistance(point: core.math.Vec2, start: core.math.Vec2, end: core.math.Vec2) SegmentDistance {
    const dx = end.x - start.x;
    const dy = end.y - start.y;
    const length_squared = dx * dx + dy * dy;
    const parameter = if (length_squared == 0)
        @as(f32, 0)
    else
        std.math.clamp(((point.x - start.x) * dx + (point.y - start.y) * dy) / length_squared, 0, 1);
    const nearest = core.math.Vec2{ .x = start.x + dx * parameter, .y = start.y + dy * parameter };
    return .{ .distance = distance(point, nearest), .parameter = parameter };
}

fn shapeAndDiscard(
    allocator: std.mem.Allocator,
    atoms: []const core.math.Vec2,
    bonds: []const LigandBond,
    pockets: []const f32,
) !void {
    const shape = try shapeAroundLigand(allocator, atoms, bonds, pockets, 0);
    allocator.free(shape);
}

/// Allocator-owned translation of upstream's `groupResiduesInSSEs` result.
/// Groups and members retain chain-lexicographic/residue-number order.
pub const SecondaryStructures = struct {
    allocator: std.mem.Allocator,
    offsets: []u32,
    residues: []core.ids.ResidueId,

    pub fn deinit(self: *SecondaryStructures) void {
        self.allocator.free(self.residues);
        self.allocator.free(self.offsets);
        self.* = undefined;
    }

    pub fn count(self: SecondaryStructures) usize {
        return self.offsets.len - 1;
    }

    pub fn members(self: SecondaryStructures, index: usize) []const core.ids.ResidueId {
        std.debug.assert(index < self.count());
        return self.residues[self.offsets[index]..self.offsets[index + 1]];
    }
};

const ResidueOrderContext = struct {
    residues: []const model.Residue,
    chains: []const u8,
};

/// Divide residues into secondary-structure runs. Upstream first orders chains
/// lexicographically and each chain by residue number, then starts a new run
/// after a gap greater than three or for blank/empty chain identifiers.
pub fn groupSecondaryStructures(
    allocator: std.mem.Allocator,
    residues: []const model.Residue,
    chain_bytes: []const u8,
) core.errors.Error!SecondaryStructures {
    const ordered = allocator.alloc(core.ids.ResidueId, residues.len) catch return error.OutOfMemory;
    errdefer allocator.free(ordered);
    for (ordered, 0..) |*id, index| id.* = core.ids.ResidueId.fromIndex(@intCast(index));
    std.mem.sortUnstable(
        core.ids.ResidueId,
        ordered,
        ResidueOrderContext{ .residues = residues, .chains = chain_bytes },
        residueLessThan,
    );

    const all_offsets = allocator.alloc(u32, residues.len + 1) catch return error.OutOfMemory;
    defer allocator.free(all_offsets);
    var group_count: usize = 0;
    var previous: ?model.Residue = null;
    for (ordered, 0..) |id, ordered_index| {
        const residue = residues[id.index()];
        const chain = try residueChain(residue, chain_bytes);
        const starts_group = if (previous) |prior|
            !std.mem.eql(u8, try residueChain(prior, chain_bytes), chain) or
                residue.residue_number -| prior.residue_number > 3 or
                chain.len == 0 or std.mem.eql(u8, chain, " ")
        else
            true;
        if (starts_group) {
            all_offsets[group_count] = @intCast(ordered_index);
            group_count += 1;
        }
        previous = residue;
    }
    all_offsets[group_count] = @intCast(ordered.len);
    const offsets = allocator.dupe(u32, all_offsets[0 .. group_count + 1]) catch return error.OutOfMemory;
    return .{ .allocator = allocator, .offsets = offsets, .residues = ordered };
}

fn residueLessThan(
    context: ResidueOrderContext,
    first_id: core.ids.ResidueId,
    second_id: core.ids.ResidueId,
) bool {
    const first = context.residues[first_id.index()];
    const second = context.residues[second_id.index()];
    const first_chain = residueChain(first, context.chains) catch unreachable;
    const second_chain = residueChain(second, context.chains) catch unreachable;
    return switch (std.mem.order(u8, first_chain, second_chain)) {
        .lt => true,
        .gt => false,
        .eq => if (first.residue_number != second.residue_number)
            first.residue_number < second.residue_number
        else
            first.id.index() < second.id.index(),
    };
}

fn residueChain(residue: model.Residue, chain_bytes: []const u8) core.errors.Error![]const u8 {
    const start = residue.chain_start;
    const end = @as(usize, start) + residue.chain_len;
    if (end > chain_bytes.len) return error.InvalidMapping;
    return chain_bytes[start..end];
}

fn residueDistance(
    start: f32,
    increment: f32,
    target: core.ids.ResidueId,
    secondary_structure: []const core.ids.ResidueId,
    residues: []const model.Residue,
) core.errors.Error!f32 {
    if (!std.math.isFinite(start) or !std.math.isFinite(increment)) return error.InvalidCoordinate;
    var total = start;
    var previous: ?model.Residue = null;
    for (secondary_structure) |id| {
        if (id.index() >= residues.len) return error.InvalidMapping;
        const residue = residues[id.index()];
        if (previous) |last| {
            const gap: f32 = @floatFromInt(@as(i64, residue.residue_number) - @as(i64, last.residue_number));
            total += increment * @max(@as(f32, 1), 1 + (gap - 1) * 0.8);
        }
        if (id == target) return total;
        previous = residue;
    }
    return error.InvalidMapping;
}

fn shapeIndex(shape_len: usize, position: f32) core.errors.Error!usize {
    if (shape_len == 0) return error.InvalidMapping;
    if (!std.math.isFinite(position)) return error.InvalidCoordinate;
    const normalized = position - @floor(position);
    const scaled = @as(f32, @floatFromInt(shape_len)) * normalized;
    return @min(@as(usize, @intFromFloat(scaled)), shape_len - 1);
}

const CrownPlacementContext = struct {
    allocator: std.mem.Allocator,
    atoms: []model.Atom,
    residues: []const model.Residue,
    interactions: []const model.ResidueInteraction,
    ligand_atoms: []const core.ids.AtomId,
    coordinates_set: []bool,
    shape: []const core.math.Vec2,
    shape_number: u32,
    penalties: []bool,

    fn validate(self: CrownPlacementContext) core.errors.Error!void {
        if (self.atoms.len != self.coordinates_set.len or self.shape.len != self.penalties.len or self.shape.len == 0) {
            return error.InvalidMapping;
        }
        for (self.ligand_atoms) |atom| if (atom.index() >= self.atoms.len) return error.InvalidMapping;
    }
};

/// Place one ordered secondary structure on a crown using upstream's exhaustive
/// start/direction scoring and shared occupancy mask. Residues that collide
/// with an earlier structure remain unset for placement on the next crown.
fn placeSecondaryStructure(
    secondary_structure: []const core.ids.ResidueId,
    context: CrownPlacementContext,
) core.errors.Error!usize {
    try context.validate();
    var to_place: usize = 0;
    for (secondary_structure) |id| {
        if (id.index() >= context.residues.len) return error.InvalidMapping;
        const atom = context.residues[id.index()].atom;
        if (atom.index() >= context.atoms.len) return error.InvalidMapping;
        to_place += @intFromBool(!context.coordinates_set[atom.index()]);
    }
    if (to_place == 0) return 0;

    const increment_distance = 5 / @as(f32, @floatFromInt(context.shape.len));
    var best_start: f32 = 0;
    var best_increment: f32 = -increment_distance;
    var best_score = std.math.inf(f32);
    var start: f32 = 0;
    while (start < 1) : (start += 0.004) {
        var increment = -increment_distance;
        while (increment <= increment_distance) : (increment += increment_distance) {
            if (increment == 0) continue;
            const score = try scoreSecondaryStructure(secondary_structure, start, increment, context);
            if (score < best_score) {
                best_score = score;
                best_start = start;
                best_increment = increment;
            }
        }
    }

    const outliers = context.allocator.alloc(bool, secondary_structure.len) catch return error.OutOfMemory;
    defer context.allocator.free(outliers);
    var placed: usize = 0;
    for (secondary_structure, outliers) |id, *outlier| {
        const residue = context.residues[id.index()];
        if (context.coordinates_set[residue.atom.index()]) {
            outlier.* = false;
            continue;
        }
        const position = try residueDistance(best_start, best_increment, id, secondary_structure, context.residues);
        const index = try shapeIndex(context.shape.len, position);
        outlier.* = context.penalties[index];
        if (outlier.*) continue;
        context.atoms[residue.atom.index()].coordinates = context.shape[index];
        placed += 1;
    }

    try markSecondaryStructure(secondary_structure, outliers, best_start, best_increment, context);
    for (secondary_structure, outliers) |id, outlier| {
        const residue = context.residues[id.index()];
        if (!context.coordinates_set[residue.atom.index()] and !outlier) context.coordinates_set[residue.atom.index()] = true;
    }
    return to_place - placed;
}

fn scoreSecondaryStructure(
    secondary_structure: []const core.ids.ResidueId,
    start: f32,
    increment: f32,
    context: CrownPlacementContext,
) core.errors.Error!f32 {
    var score: f32 = 0;
    var previous_state: ?i2 = null;
    var previous_coordinates: core.math.Vec2 = .{};
    for (secondary_structure) |id| {
        const residue = context.residues[id.index()];
        const atom_index = residue.atom.index();
        const position = try residueDistance(start, increment, id, secondary_structure, context.residues);
        const index = try shapeIndex(context.shape.len, position);
        const state: i2 = if (context.coordinates_set[atom_index]) -1 else if (context.penalties[index]) 1 else 0;
        const coordinates = if (state == -1) context.atoms[atom_index].coordinates else context.shape[index];
        if (state != -1) score += try scoreResiduePosition(residue, coordinates, context);
        if (previous_state) |prior| {
            if (state != prior) score += distanceSquared(coordinates, previous_coordinates) * 400;
        }
        previous_state = state;
        previous_coordinates = coordinates;
    }
    return score;
}

fn scoreResiduePosition(
    residue: model.Residue,
    position: core.math.Vec2,
    context: CrownPlacementContext,
) core.errors.Error!f32 {
    var score: f32 = 0;
    var target_count: usize = 0;
    for (context.interactions) |interaction| {
        const target = if (interaction.start == residue.atom)
            interaction.end
        else if (interaction.end == residue.atom)
            interaction.start
        else
            continue;
        if (target.index() >= context.atoms.len) return error.InvalidMapping;
        if (!context.coordinates_set[target.index()]) continue;
        score += scoreTarget(position, target, 1, context);
        target_count += 1;
    }
    if (target_count == 0 and residue.closest_ligand_atom.isValid()) {
        if (residue.closest_ligand_atom.index() >= context.atoms.len) return error.InvalidMapping;
        score += scoreTarget(position, residue.closest_ligand_atom, 0.2, context);
    }
    return score;
}

fn scoreTarget(
    position: core.math.Vec2,
    target: core.ids.AtomId,
    interaction_scaling: f32,
    context: CrownPlacementContext,
) f32 {
    const target_coordinates = context.atoms[target.index()].coordinates;
    var clashing: f32 = 0;
    for (context.ligand_atoms) |ligand| {
        if (ligand == target) continue;
        if (pointSegmentDistance(context.atoms[ligand.index()].coordinates, position, target_coordinates).distance < 40) {
            clashing += 1;
        }
    }
    const expected_distance = @as(f32, @floatFromInt(context.shape_number)) * 50;
    return interaction_scaling * (0.01 * (distanceSquared(target_coordinates, position) - expected_distance * expected_distance) +
        clashing * 100);
}

fn markSecondaryStructure(
    secondary_structure: []const core.ids.ResidueId,
    outliers: []const bool,
    start: f32,
    increment: f32,
    context: CrownPlacementContext,
) core.errors.Error!void {
    if (secondary_structure.len != outliers.len) return error.InvalidMapping;
    const padding = @abs(increment) * 0.5;
    var previous_position: ?f32 = null;
    for (secondary_structure, outliers) |id, outlier| {
        const residue = context.residues[id.index()];
        if (context.coordinates_set[residue.atom.index()] or outlier) {
            previous_position = null;
            continue;
        }
        const position = try residueDistance(start, increment, id, secondary_structure, context.residues);
        markShapeRange(context.penalties, try shapeIndex(context.shape.len, position - padding), try shapeIndex(context.shape.len, position + padding));
        if (previous_position) |prior| {
            const range_start = if (increment < 0) position else prior;
            const range_end = if (increment < 0) prior else position;
            markShapeRange(context.penalties, try shapeIndex(context.shape.len, range_start), try shapeIndex(context.shape.len, range_end));
        }
        previous_position = position;
    }
}

fn markShapeRange(penalties: []bool, start: usize, end: usize) void {
    var index = start;
    while (index != end) : (index = (index + 1) % penalties.len) penalties[index] = true;
}

fn testAtom(index: u32, coordinates: core.math.Vec2) model.Atom {
    return .{
        .id = core.ids.AtomId.fromIndex(index),
        .input_index = index,
        .atomic_number = .carbon,
        .coordinates = coordinates,
    };
}

fn placeAroundAndDiscard(
    allocator: std.mem.Allocator,
    source_atoms: []const model.Atom,
    bonds: []const model.Bond,
    residue_input: []const model.Residue,
    chains: []const u8,
    interactions: []const model.ResidueInteraction,
) !void {
    const atoms = try allocator.dupe(model.Atom, source_atoms);
    defer allocator.free(atoms);
    try placeAroundLigand(allocator, atoms, bonds, residue_input, chains, interactions);
}

test "secondary structures follow chain order, residue gaps, and blank-chain splits" {
    const chains = "BBAAA  ";
    const input = [_]model.Residue{
        .{ .id = core.ids.ResidueId.fromIndex(0), .atom = core.ids.AtomId.fromIndex(0), .chain_start = 0, .chain_len = 1, .residue_number = 9 },
        .{ .id = core.ids.ResidueId.fromIndex(1), .atom = core.ids.AtomId.fromIndex(1), .chain_start = 1, .chain_len = 1, .residue_number = 2 },
        .{ .id = core.ids.ResidueId.fromIndex(2), .atom = core.ids.AtomId.fromIndex(2), .chain_start = 2, .chain_len = 1, .residue_number = 1 },
        .{ .id = core.ids.ResidueId.fromIndex(3), .atom = core.ids.AtomId.fromIndex(3), .chain_start = 3, .chain_len = 1, .residue_number = 3 },
        .{ .id = core.ids.ResidueId.fromIndex(4), .atom = core.ids.AtomId.fromIndex(4), .chain_start = 4, .chain_len = 1, .residue_number = 8 },
        .{ .id = core.ids.ResidueId.fromIndex(5), .atom = core.ids.AtomId.fromIndex(5), .chain_start = 5, .chain_len = 1, .residue_number = 1 },
        .{ .id = core.ids.ResidueId.fromIndex(6), .atom = core.ids.AtomId.fromIndex(6), .chain_start = 6, .chain_len = 1, .residue_number = 2 },
    };
    var groups = try groupSecondaryStructures(std.testing.allocator, &input, chains);
    defer groups.deinit();
    try std.testing.expectEqual(@as(usize, 6), groups.count());
    try std.testing.expectEqualSlices(core.ids.ResidueId, &.{core.ids.ResidueId.fromIndex(5)}, groups.members(0));
    try std.testing.expectEqualSlices(core.ids.ResidueId, &.{core.ids.ResidueId.fromIndex(6)}, groups.members(1));
    try std.testing.expectEqualSlices(core.ids.ResidueId, &.{
        core.ids.ResidueId.fromIndex(2),
        core.ids.ResidueId.fromIndex(3),
    }, groups.members(2));
    try std.testing.expectEqualSlices(core.ids.ResidueId, &.{core.ids.ResidueId.fromIndex(4)}, groups.members(3));
    try std.testing.expectEqualSlices(core.ids.ResidueId, &.{core.ids.ResidueId.fromIndex(1)}, groups.members(4));
    try std.testing.expectEqualSlices(core.ids.ResidueId, &.{core.ids.ResidueId.fromIndex(0)}, groups.members(5));
}

test "ligand crowns are deterministic ordered contours at successive distances" {
    const atoms = [_]core.math.Vec2{.{}};
    const pockets = [_]f32{0};
    const first = try shapeAroundLigand(std.testing.allocator, &atoms, &.{}, &pockets, 0);
    defer std.testing.allocator.free(first);
    const repeated = try shapeAroundLigand(std.testing.allocator, &atoms, &.{}, &pockets, 0);
    defer std.testing.allocator.free(repeated);
    const second = try shapeAroundLigand(std.testing.allocator, &atoms, &.{}, &pockets, 1);
    defer std.testing.allocator.free(second);

    try std.testing.expect(first.len > 8);
    try std.testing.expectEqualSlices(core.math.Vec2, first, repeated);
    for (first) |point| try std.testing.expectApproxEqAbs(@as(f32, 60), distance(.{}, point), 1);
    for (second) |point| try std.testing.expectApproxEqAbs(@as(f32, 120), distance(.{}, point), 1);
}

test "ligand crowns include bonds and clean allocation failures" {
    const atoms = [_]core.math.Vec2{ .{ .x = -50 }, .{ .x = 50 } };
    const bonds = [_]LigandBond{.{ .start = 0, .end = 1 }};
    const pockets = [_]f32{ 0, 0 };
    const shape = try shapeAroundLigand(std.testing.allocator, &atoms, &bonds, &pockets, 0);
    defer std.testing.allocator.free(shape);
    try std.testing.expect(shape.len > 8);
    try std.testing.checkAllAllocationFailures(std.testing.allocator, shapeAndDiscard, .{ &atoms, &bonds, &pockets });
    try std.testing.expectError(error.InvalidMapping, shapeAroundLigand(std.testing.allocator, &atoms, &.{.{ .start = 0, .end = 2 }}, &pockets, 0));
}

test "secondary-structure crown spacing contracts missing residue gaps" {
    const input = [_]model.Residue{
        .{ .id = core.ids.ResidueId.fromIndex(0), .atom = core.ids.AtomId.fromIndex(0), .chain_start = 0, .chain_len = 1, .residue_number = 2 },
        .{ .id = core.ids.ResidueId.fromIndex(1), .atom = core.ids.AtomId.fromIndex(1), .chain_start = 0, .chain_len = 1, .residue_number = 5 },
        .{ .id = core.ids.ResidueId.fromIndex(2), .atom = core.ids.AtomId.fromIndex(2), .chain_start = 0, .chain_len = 1, .residue_number = 5 },
    };
    const group = [_]core.ids.ResidueId{
        core.ids.ResidueId.fromIndex(0),
        core.ids.ResidueId.fromIndex(1),
        core.ids.ResidueId.fromIndex(2),
    };
    try std.testing.expectApproxEqAbs(@as(f32, 0.1), try residueDistance(0.1, 0.2, group[0], &group, &input), 0.0001);
    try std.testing.expectApproxEqAbs(@as(f32, 0.62), try residueDistance(0.1, 0.2, group[1], &group, &input), 0.0001);
    try std.testing.expectApproxEqAbs(@as(f32, 0.82), try residueDistance(0.1, 0.2, group[2], &group, &input), 0.0001);
    try std.testing.expectEqual(@as(usize, 9), try shapeIndex(10, -0.01));
    try std.testing.expectEqual(@as(usize, 0), try shapeIndex(10, 1));
}

test "secondary structures choose interaction-facing crown positions and reserve them" {
    var atoms = [_]model.Atom{
        testAtom(0, .{ .x = -20 }),
        testAtom(1, .{ .x = 20 }),
        testAtom(2, .{}),
        testAtom(3, .{}),
    };
    const residue_input = [_]model.Residue{
        .{ .id = core.ids.ResidueId.fromIndex(0), .atom = core.ids.AtomId.fromIndex(2), .chain_start = 0, .chain_len = 1, .residue_number = 1 },
        .{ .id = core.ids.ResidueId.fromIndex(1), .atom = core.ids.AtomId.fromIndex(3), .chain_start = 0, .chain_len = 1, .residue_number = 2 },
    };
    const interactions = [_]model.ResidueInteraction{
        .{ .start = core.ids.AtomId.fromIndex(2), .end = core.ids.AtomId.fromIndex(0) },
        .{ .start = core.ids.AtomId.fromIndex(3), .end = core.ids.AtomId.fromIndex(1) },
    };
    const shape = [_]core.math.Vec2{
        .{ .x = -60 },          .{ .x = -50, .y = 30 }, .{ .x = -30, .y = 50 },  .{ .y = 60 },
        .{ .x = 30, .y = 50 },  .{ .x = 50, .y = 30 },  .{ .x = 60 },            .{ .x = 50, .y = -30 },
        .{ .x = 30, .y = -50 }, .{ .y = -60 },          .{ .x = -30, .y = -50 }, .{ .x = -50, .y = -30 },
    };
    var coordinates_set = [_]bool{ true, true, false, false };
    var penalties: [shape.len]bool = undefined;
    @memset(&penalties, false);
    const group = [_]core.ids.ResidueId{ core.ids.ResidueId.fromIndex(0), core.ids.ResidueId.fromIndex(1) };
    const remaining = try placeSecondaryStructure(&group, .{
        .allocator = std.testing.allocator,
        .atoms = &atoms,
        .residues = &residue_input,
        .interactions = &interactions,
        .ligand_atoms = &.{ core.ids.AtomId.fromIndex(0), core.ids.AtomId.fromIndex(1) },
        .coordinates_set = &coordinates_set,
        .shape = &shape,
        .shape_number = 1,
        .penalties = &penalties,
    });
    try std.testing.expectEqual(@as(usize, 0), remaining);
    try std.testing.expect(coordinates_set[2] and coordinates_set[3]);
    try std.testing.expect(atoms[2].coordinates.x != atoms[3].coordinates.x or atoms[2].coordinates.y != atoms[3].coordinates.y);
    try std.testing.expect(std.mem.indexOfScalar(bool, &penalties, true) != null);
}

test "residue crown placement is deterministic and allocation-clean" {
    const source_atoms = [_]model.Atom{
        testAtom(0, .{ .x = -25 }),
        testAtom(1, .{ .x = 25 }),
        testAtom(2, .{}),
    };
    const bonds = [_]model.Bond{.{
        .id = core.ids.BondId.fromIndex(0),
        .input_index = 0,
        .start = core.ids.AtomId.fromIndex(0),
        .end = core.ids.AtomId.fromIndex(1),
        .input_order = .single,
        .effective_order = .single,
    }};
    const residue_input = [_]model.Residue{.{
        .id = core.ids.ResidueId.fromIndex(0),
        .atom = core.ids.AtomId.fromIndex(2),
        .chain_start = 0,
        .chain_len = 1,
        .residue_number = 7,
        .closest_ligand_atom = core.ids.AtomId.fromIndex(1),
    }};
    const interactions = [_]model.ResidueInteraction{.{
        .start = core.ids.AtomId.fromIndex(2),
        .end = core.ids.AtomId.fromIndex(1),
    }};
    var first = source_atoms;
    var second = source_atoms;
    try placeAroundLigand(std.testing.allocator, &first, &bonds, &residue_input, "A", &interactions);
    try placeAroundLigand(std.testing.allocator, &second, &bonds, &residue_input, "A", &interactions);
    try std.testing.expectEqual(first[2].coordinates, second[2].coordinates);
    try std.testing.expect(distance(first[2].coordinates, first[1].coordinates) > 40);
    try std.testing.checkAllAllocationFailures(std.testing.allocator, placeAroundAndDiscard, .{
        &source_atoms,
        &bonds,
        &residue_input,
        "A",
        &interactions,
    });
}
