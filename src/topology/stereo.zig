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
