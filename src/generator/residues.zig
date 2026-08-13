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
