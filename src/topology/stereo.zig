const std = @import("std");
const core = @import("core");
const model = @import("model");

pub const Priority = enum { left, right, tied };

/// Compare two substituent branches using the CIP breadth-layer rule. Atomic
/// numbers are compared first; subsequent shells are sorted descending, with
/// duplicate ghost entries for multiple bonds and no repeated traversal of a
/// previously visited atom. This is the value-owned counterpart of CIPAtom.
pub fn compareBranches(
    allocator: std.mem.Allocator,
    atoms: []const model.Atom,
    bonds: []const model.Bond,
    graph: anytype,
    center: core.ids.AtomId,
    left: core.ids.AtomId,
    right: core.ids.AtomId,
) core.errors.Error!Priority {
    const left_number = @backingInt(atoms[left.index()].atomic_number);
    const right_number = @backingInt(atoms[right.index()].atomic_number);
    if (left_number > right_number) return .left;
    if (right_number > left_number) return .right;

    var left_branch = try Branch.init(allocator, atoms.len, center, left);
    defer left_branch.deinit();
    var right_branch = try Branch.init(allocator, atoms.len, center, right);
    defer right_branch.deinit();
    while (left_branch.frontier.items.len != 0 or right_branch.frontier.items.len != 0) {
        const left_shell = try left_branch.expand(allocator, atoms, bonds, graph);
        defer allocator.free(left_shell);
        const right_shell = try right_branch.expand(allocator, atoms, bonds, graph);
        defer allocator.free(right_shell);
        const comparison = compareShells(left_shell, right_shell);
        if (comparison != .tied) return comparison;
    }
    return .tied;
}

const Branch = struct {
    const Node = struct { atom: core.ids.AtomId, parent: core.ids.AtomId };

    allocator: std.mem.Allocator,
    visited: []bool,
    frontier: std.ArrayList(Node),

    fn init(allocator: std.mem.Allocator, count: usize, center: core.ids.AtomId, first: core.ids.AtomId) core.errors.Error!Branch {
        const visited = allocator.alloc(bool, count) catch return error.OutOfMemory;
        errdefer allocator.free(visited);
        @memset(visited, false);
        visited[center.index()] = true;
        visited[first.index()] = true;
        var frontier: std.ArrayList(Node) = .empty;
        frontier.append(allocator, .{ .atom = first, .parent = center }) catch return error.OutOfMemory;
        return .{ .allocator = allocator, .visited = visited, .frontier = frontier };
    }

    fn deinit(self: *Branch) void {
        self.frontier.deinit(self.allocator);
        self.allocator.free(self.visited);
    }

    fn expand(self: *Branch, allocator: std.mem.Allocator, atoms: []const model.Atom, bonds: []const model.Bond, graph: anytype) core.errors.Error![]u32 {
        var shell: std.ArrayList(u32) = .empty;
        defer shell.deinit(allocator);
        var next: std.ArrayList(Node) = .empty;
        errdefer next.deinit(allocator);
        for (self.frontier.items) |node| {
            for (graph.neighbors(node.atom), graph.incidentBonds(node.atom)) |neighbor, bond_id| {
                const bond = bonds[bond_id.index()];
                const number = @backingInt(atoms[neighbor.index()].atomic_number);
                const duplicate_count: u32 = switch (bond.effective_order) {
                    .double => 1,
                    .triple => 2,
                    else => 0,
                };
                for (0..duplicate_count) |_| shell.append(allocator, number) catch return error.OutOfMemory;
                if (neighbor == node.parent) continue;
                if (self.visited[neighbor.index()]) {
                    shell.append(allocator, number) catch return error.OutOfMemory;
                    continue;
                }
                self.visited[neighbor.index()] = true;
                shell.append(allocator, number) catch return error.OutOfMemory;
                next.append(allocator, .{ .atom = neighbor, .parent = node.atom }) catch return error.OutOfMemory;
            }
        }
        std.mem.sort(u32, shell.items, {}, std.sort.desc(u32));
        self.frontier.deinit(self.allocator);
        self.frontier = next;
        return shell.toOwnedSlice(allocator) catch return error.OutOfMemory;
    }
};

fn compareShells(left: []const u32, right: []const u32) Priority {
    const shared = @min(left.len, right.len);
    for (left[0..shared], right[0..shared]) |left_number, right_number| {
        if (left_number > right_number) return .left;
        if (right_number > left_number) return .right;
    }
    if (left.len > right.len) return .left;
    if (right.len > left.len) return .right;
    return .tied;
}

pub fn firstNeighbor(
    allocator: std.mem.Allocator,
    atoms: []const model.Atom,
    bonds: []const model.Bond,
    graph: anytype,
    endpoint: core.ids.AtomId,
    other_endpoint: core.ids.AtomId,
) core.errors.Error!?core.ids.AtomId {
    var best: ?core.ids.AtomId = null;
    for (graph.neighbors(endpoint)) |neighbor| {
        if (neighbor == other_endpoint) continue;
        if (best) |current| {
            if (try compareBranches(allocator, atoms, bonds, graph, endpoint, neighbor, current) == .left) best = neighbor;
        } else best = neighbor;
    }
    return best;
}

pub fn absoluteBondStereo(
    allocator: std.mem.Allocator,
    atoms: []const model.Atom,
    bonds: []const model.Bond,
    graph: anytype,
    membership: anytype,
    bond_id: core.ids.BondId,
) core.errors.Error!core.chemistry.BondStereo {
    const bond = bonds[bond_id.index()];
    if (!isStereoBond(membership, bond)) return .unspecified;
    const start_neighbor = (try firstNeighbor(allocator, atoms, bonds, graph, bond.start, bond.end)) orelse return .unspecified;
    const end_neighbor = (try firstNeighbor(allocator, atoms, bonds, graph, bond.end, bond.start)) orelse return .unspecified;
    if (!bond.stereo_atom_a.isValid() or !bond.stereo_atom_b.isValid()) return .unspecified;
    var invert = bond.stereo_atom_a != start_neighbor and bond.stereo_atom_a != end_neighbor;
    if (bond.stereo_atom_b != start_neighbor and bond.stereo_atom_b != end_neighbor) invert = !invert;
    var is_z = bond.stereo == .cis or bond.stereo == .z;
    if (invert) is_z = !is_z;
    return if (is_z) .z else .e;
}

pub fn isStereoBond(membership: anytype, bond: model.Bond) bool {
    if (bond.effective_order != .double or bond.stereo == .unspecified) return false;
    for (membership.bondRings(bond.id)) |ring| {
        if (membership.atoms(ring).len < 9) return false;
    }
    return true;
}

pub fn geometryMatches(
    atoms: []const model.Atom,
    bond: model.Bond,
    start_neighbor: core.ids.AtomId,
    end_neighbor: core.ids.AtomId,
    absolute: core.chemistry.BondStereo,
) bool {
    if (absolute != .z and absolute != .e) return true;
    const start = atoms[bond.start.index()].coordinates;
    const end = atoms[bond.end.index()].coordinates;
    const a = atoms[start_neighbor.index()].coordinates;
    const b = atoms[end_neighbor.index()].coordinates;
    const side_a = (end.x - start.x) * (a.y - start.y) - (end.y - start.y) * (a.x - start.x);
    const side_b = (end.x - start.x) * (b.y - start.y) - (end.y - start.y) * (b.x - start.x);
    const same_side = side_a * side_b > 0;
    return same_side == (absolute == .z);
}

const Substituent = union(enum) {
    atom: core.ids.AtomId,
    implicit_hydrogen,
};

const RankedNeighbor = struct {
    atom: core.ids.AtomId,
    priority: u2,
};

/// Convert caller-relative tetrahedral chirality to an absolute CIP label.
/// The relative order is atom-a, atom-b, the remaining substituent, then the
/// atom being looked from. An omitted substituent is the one implicit H.
pub fn absoluteAtomStereo(
    allocator: std.mem.Allocator,
    atoms: []const model.Atom,
    bonds: []const model.Bond,
    graph: anytype,
    center: core.ids.AtomId,
) core.errors.Error!core.chemistry.AtomStereo {
    if (center.index() >= atoms.len) return error.InvalidAtomIndex;
    const atom = atoms[center.index()];
    if (atom.stereo == .unspecified or atom.stereo == .r or atom.stereo == .s) return atom.stereo;
    const degree = graph.degree(center);
    if (degree != 3 and degree != 4) return error.InvalidStereo;
    const references = [_]core.ids.AtomId{ atom.stereo_atom_a, atom.stereo_atom_b, atom.stereo_looking_from };
    for (references, 0..) |reference, index| {
        if (!reference.isValid() or reference.index() >= atoms.len or reference == center or
            !containsAtom(graph.neighbors(center), reference)) return error.InvalidStereo;
        for (references[0..index]) |previous| if (reference == previous) return error.InvalidStereo;
    }

    var ordered: [4]Substituent = undefined;
    ordered[0] = .{ .atom = atom.stereo_atom_a };
    ordered[1] = .{ .atom = atom.stereo_atom_b };
    ordered[3] = .{ .atom = atom.stereo_looking_from };
    var remaining: ?core.ids.AtomId = null;
    for (graph.neighbors(center)) |neighbor| {
        if (containsAtom(&references, neighbor)) continue;
        if (remaining != null) return error.InvalidStereo;
        remaining = neighbor;
    }
    ordered[2] = if (remaining) |neighbor| .{ .atom = neighbor } else .implicit_hydrogen;
    if (degree == 4 and remaining == null) return error.InvalidStereo;

    var priorities: [4]u2 = undefined;
    for (0..4) |i| {
        var rank: u2 = 0;
        for (0..4) |j| {
            if (i == j) continue;
            switch (try compareSubstituents(allocator, atoms, bonds, graph, center, ordered[j], ordered[i])) {
                .left => rank += 1,
                .right => {},
                .tied => return error.InvalidStereo,
            }
        }
        priorities[i] = rank;
    }
    var inversions: u32 = 0;
    for (priorities, 0..) |left, i| for (priorities[i + 1 ..]) |right| {
        inversions += @intFromBool(left > right);
    };
    var is_r = atom.stereo == .counter_clockwise;
    if (inversions % 2 == 1) is_r = !is_r;
    return if (is_r) .r else .s;
}

fn compareSubstituents(
    allocator: std.mem.Allocator,
    atoms: []const model.Atom,
    bonds: []const model.Bond,
    graph: anytype,
    center: core.ids.AtomId,
    left: Substituent,
    right: Substituent,
) core.errors.Error!Priority {
    if (left == .implicit_hydrogen) return if (right == .implicit_hydrogen) .tied else .right;
    if (right == .implicit_hydrogen) return .left;
    return compareBranches(allocator, atoms, bonds, graph, center, left.atom, right.atom);
}

fn containsAtom(haystack: []const core.ids.AtomId, needle: core.ids.AtomId) bool {
    return std.mem.indexOfScalar(core.ids.AtomId, haystack, needle) != null;
}

/// Assign a deterministic wedge/hash pair which realizes each absolute atom
/// descriptor in the current 2D geometry. Four-coordinate centers receive a
/// solid/hash pair; three-coordinate centers receive one display bond and use
/// the implicit substituent on the opposite side of the drawing plane.
pub fn writeAtomBondDisplays(
    allocator: std.mem.Allocator,
    atoms: []model.Atom,
    bonds: []model.Bond,
    graph: anytype,
    membership: anytype,
) core.errors.Error!void {
    for (bonds) |*bond| bond.display = .none;
    for (atoms) |*atom| {
        const absolute = try absoluteAtomStereo(allocator, atoms, bonds, graph, atom.id);
        atom.stereo = absolute;
        if (absolute != .r and absolute != .s) continue;
        const neighbors = graph.neighbors(atom.id);
        if (neighbors.len != 3 and neighbors.len != 4) return error.InvalidStereo;

        var ranked: [4]RankedNeighbor = undefined;
        var display_order: [4]RankedNeighbor = undefined;
        var display_weights: [4]f32 = undefined;
        for (neighbors, 0..) |neighbor, i| ranked[i].atom = neighbor;
        for (neighbors, 0..) |_, i| {
            var rank: u2 = 0;
            for (neighbors, 0..) |_, j| {
                if (i == j) continue;
                switch (try compareBranches(allocator, atoms, bonds, graph, atom.id, ranked[j].atom, ranked[i].atom)) {
                    .left => rank += 1,
                    .right => {},
                    .tied => return error.InvalidStereo,
                }
            }
            if (neighbors.len == 3) rank += 0;
            ranked[i].priority = rank;
        }
        std.mem.sort(@TypeOf(ranked[0]), ranked[0..neighbors.len], {}, priorityLessThan);

        for (neighbors, 0..) |neighbor, i| {
            display_order[i] = .{ .atom = neighbor, .priority = 0 };
            const branch = try graph.reachableExcluding(allocator, neighbor, atom.id);
            defer allocator.free(branch);
            const bond_id = bondBetween(graph, atom.id, neighbor) orelse return error.InvalidStereo;
            const bond = bonds[bond_id.index()];
            var weight: f32 = @floatFromInt(branch.len);
            if (bond.effective_order == .double) weight -= 0.25;
            if (atom.atomic_number == .sulfur and bond.effective_order == .double) weight += 2000;
            if (membership.bondRings(bond_id).len != 0) weight += 500;
            if (atoms[neighbor.index()].atomic_number == .carbon) weight += 0.5;
            if (atoms[neighbor.index()].atomic_number == .hydrogen) weight -= 0.5;
            if (atoms[neighbor.index()].stereo != .unspecified) weight += 10000;
            if (atom.cross_layout and graph.degree(neighbor) > 1) weight += 200;
            for (graph.incidentBonds(neighbor)) |neighbor_bond| if (bonds[neighbor_bond.index()].effective_order == .double) {
                weight += 100;
                break;
            };
            display_weights[i] = weight;
        }
        sortByDisplayWeight(display_order[0..neighbors.len], display_weights[0..neighbors.len]);

        const first_bond = bondBetween(graph, atom.id, display_order[0].atom) orelse return error.InvalidStereo;
        bonds[first_bond.index()].display = displayFrom(atom.id, bonds[first_bond.index()], true);
        if (neighbors.len == 4) {
            const second_bond = bondBetween(graph, atom.id, display_order[1].atom) orelse return error.InvalidStereo;
            bonds[second_bond.index()].display = displayFrom(atom.id, bonds[second_bond.index()], false);
        }
        if (displayChirality(atoms, bonds, graph, atom.id, ranked[0..neighbors.len]) != absolute) {
            bonds[first_bond.index()].display = invertDisplay(bonds[first_bond.index()].display);
            if (neighbors.len == 4) {
                const second_bond = bondBetween(graph, atom.id, display_order[1].atom).?;
                bonds[second_bond.index()].display = invertDisplay(bonds[second_bond.index()].display);
            }
        }
    }
}

fn priorityLessThan(_: void, left: RankedNeighbor, right: RankedNeighbor) bool {
    return left.priority < right.priority;
}

fn sortByDisplayWeight(neighbors: []RankedNeighbor, weights: []f32) void {
    for (1..neighbors.len) |index| {
        var cursor = index;
        while (cursor > 0 and weights[cursor] < weights[cursor - 1]) : (cursor -= 1) {
            std.mem.swap(f32, &weights[cursor], &weights[cursor - 1]);
            std.mem.swap(RankedNeighbor, &neighbors[cursor], &neighbors[cursor - 1]);
        }
    }
}

fn bondBetween(graph: anytype, center: core.ids.AtomId, neighbor: core.ids.AtomId) ?core.ids.BondId {
    for (graph.neighbors(center), graph.incidentBonds(center)) |candidate, bond| if (candidate == neighbor) return bond;
    return null;
}

fn displayFrom(center: core.ids.AtomId, bond: model.Bond, solid: bool) core.chemistry.BondDisplay {
    const forward = bond.start == center;
    return if (solid)
        (if (forward) .solid_forward else .solid_reverse)
    else
        (if (forward) .hashed_forward else .hashed_reverse);
}

fn invertDisplay(display: core.chemistry.BondDisplay) core.chemistry.BondDisplay {
    return switch (display) {
        .solid_forward => .hashed_forward,
        .solid_reverse => .hashed_reverse,
        .hashed_forward => .solid_forward,
        .hashed_reverse => .solid_reverse,
        .none => .none,
    };
}

fn displayChirality(atoms: []const model.Atom, bonds: []const model.Bond, graph: anytype, center: core.ids.AtomId, ranked: anytype) core.chemistry.AtomStereo {
    var vectors: [4]struct { x: f32, y: f32, z: f32 } = undefined;
    const origin = atoms[center.index()].coordinates;
    for (ranked, 0..) |entry, i| {
        const position = atoms[entry.atom.index()].coordinates;
        const bond = bonds[bondBetween(graph, center, entry.atom).?.index()];
        vectors[i] = .{ .x = position.x - origin.x, .y = position.y - origin.y, .z = displayZ(center, bond) };
    }
    if (ranked.len == 3) {
        vectors[3] = .{ .x = 0, .y = 0, .z = -(vectors[0].z + vectors[1].z + vectors[2].z) };
    }
    const a = subtract3(vectors[0], vectors[3]);
    const b = subtract3(vectors[1], vectors[3]);
    const c = subtract3(vectors[2], vectors[3]);
    const determinant = a.x * (b.y * c.z - b.z * c.y) -
        a.y * (b.x * c.z - b.z * c.x) + a.z * (b.x * c.y - b.y * c.x);
    return if (determinant < 0) .r else .s;
}

fn subtract3(a: anytype, b: @TypeOf(a)) @TypeOf(a) {
    return .{ .x = a.x - b.x, .y = a.y - b.y, .z = a.z - b.z };
}

fn displayZ(center: core.ids.AtomId, bond: model.Bond) f32 {
    return switch (bond.display) {
        .solid_forward => if (bond.start == center) 1 else -1,
        .solid_reverse => if (bond.end == center) 1 else -1,
        .hashed_forward => if (bond.start == center) -1 else 1,
        .hashed_reverse => if (bond.end == center) -1 else 1,
        .none => 0,
    };
}

/// Propagate display-bond depth over each component, matching assignPseudoZ.
pub fn assignPseudoZ(allocator: std.mem.Allocator, atoms: []const model.Atom, bonds: []const model.Bond, graph: anytype) core.errors.Error![]f32 {
    const result = allocator.alloc(f32, atoms.len) catch return error.OutOfMemory;
    errdefer allocator.free(result);
    @memset(result, 0);
    const visited = allocator.alloc(bool, atoms.len) catch return error.OutOfMemory;
    defer allocator.free(visited);
    @memset(visited, false);
    const queue = allocator.alloc(core.ids.AtomId, atoms.len) catch return error.OutOfMemory;
    defer allocator.free(queue);
    for (atoms) |atom| {
        if (atom.hidden or visited[atom.id.index()]) continue;
        queue[0] = atom.id;
        visited[atom.id.index()] = true;
        var head: usize = 0;
        var tail: usize = 1;
        while (head < tail) : (head += 1) {
            const current = queue[head];
            for (graph.neighbors(current), graph.incidentBonds(current)) |neighbor, bond_id| {
                if (visited[neighbor.index()]) continue;
                visited[neighbor.index()] = true;
                result[neighbor.index()] = result[current.index()] + displayZ(current, bonds[bond_id.index()]);
                queue[tail] = neighbor;
                tail += 1;
            }
        }
    }
    return result;
}

/// Orientation preference used when choosing between a component and its
/// horizontal mirror. Each peptide alpha carbon contributes upstream's fixed
/// ±100 score according to the N-to-carbonyl direction.
pub fn peptideFlipScore(atoms: []const model.Atom, bonds: []const model.Bond, graph: anytype) f32 {
    var score: f32 = 0;
    for (atoms) |alpha| {
        if (alpha.hidden or alpha.atomic_number != .carbon or isCarbonylCarbon(atoms, bonds, graph, alpha.id)) continue;
        var amino: ?core.ids.AtomId = null;
        var carbonyl: ?core.ids.AtomId = null;
        for (graph.neighbors(alpha.id)) |neighbor| {
            if (atoms[neighbor.index()].atomic_number == .nitrogen) amino = neighbor;
            if (isCarbonylCarbon(atoms, bonds, graph, neighbor)) carbonyl = neighbor;
        }
        if (amino != null and carbonyl != null) {
            score += if (atoms[amino.?.index()].coordinates.x > atoms[carbonyl.?.index()].coordinates.x) -100 else 100;
        }
    }
    return score;
}

fn isCarbonylCarbon(atoms: []const model.Atom, bonds: []const model.Bond, graph: anytype, atom: core.ids.AtomId) bool {
    if (atoms[atom.index()].atomic_number != .carbon) return false;
    for (graph.neighbors(atom), graph.incidentBonds(atom)) |neighbor, bond| {
        if (atoms[neighbor.index()].atomic_number == .oxygen and bonds[bond.index()].effective_order == .double) return true;
    }
    return false;
}
