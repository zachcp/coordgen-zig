const std = @import("std");
const core = @import("core");
const model = @import("model");

pub const AtomId = core.ids.AtomId;
pub const BondId = core.ids.BondId;
pub const MoleculeId = core.ids.MoleculeId;

/// Allocator-owned graph indexing after upstream's skip/zero-order/hidden
/// filtering. Neighbor and incident-bond slices are parallel and retain bond
/// input order, matching assignBondsAndNeighbors.
pub const Graph = struct {
    allocator: std.mem.Allocator,
    adjacency_offsets: []u32,
    adjacent_atoms: []AtomId,
    incident_bonds: []BondId,
    structural_bonds: []BondId,
    component_of: []MoleculeId,
    component_offsets: []u32,
    component_atoms: []AtomId,
    component_count: u32,

    pub fn init(
        allocator: std.mem.Allocator,
        atoms: []const model.Atom,
        bonds: []const model.Bond,
    ) core.errors.Error!Graph {
        if (atoms.len >= std.math.maxInt(u32) or bonds.len >= std.math.maxInt(u32)) {
            return error.TooManyItems;
        }
        for (atoms, 0..) |atom, index| {
            if (atom.id.index() != index) return error.InvalidMapping;
        }
        for (bonds, 0..) |bond, index| {
            if (bond.id.index() != index) return error.InvalidMapping;
            const start = bond.start.index();
            const end = bond.end.index();
            if (start >= atoms.len or end >= atoms.len or start == end) {
                return error.InvalidAtomIndex;
            }
        }

        const adjacency_offsets = allocator.alloc(u32, atoms.len + 1) catch return error.OutOfMemory;
        errdefer allocator.free(adjacency_offsets);
        @memset(adjacency_offsets, 0);
        var structural_bond_count: usize = 0;
        for (bonds) |bond| {
            if (!isStructuralBond(atoms, bond)) continue;
            structural_bond_count += 1;
            adjacency_offsets[bond.start.index() + 1] += 1;
            adjacency_offsets[bond.end.index() + 1] += 1;
        }
        if (structural_bond_count > std.math.maxInt(u32) / 2) return error.TooManyItems;
        for (adjacency_offsets[1..], adjacency_offsets[0 .. adjacency_offsets.len - 1]) |*count, previous| {
            count.* += previous;
        }

        const adjacency_count = adjacency_offsets[atoms.len];
        const adjacent_atoms = allocator.alloc(AtomId, adjacency_count) catch return error.OutOfMemory;
        errdefer allocator.free(adjacent_atoms);
        const incident_bonds = allocator.alloc(BondId, adjacency_count) catch return error.OutOfMemory;
        errdefer allocator.free(incident_bonds);
        const structural_bonds = allocator.alloc(BondId, structural_bond_count) catch return error.OutOfMemory;
        errdefer allocator.free(structural_bonds);
        const cursors = allocator.dupe(u32, adjacency_offsets[0..atoms.len]) catch return error.OutOfMemory;
        defer allocator.free(cursors);
        var structural_index: usize = 0;
        for (bonds) |bond| {
            if (!isStructuralBond(atoms, bond)) continue;
            structural_bonds[structural_index] = bond.id;
            structural_index += 1;
            appendEdge(cursors, adjacent_atoms, incident_bonds, bond.start, bond.end, bond.id);
            appendEdge(cursors, adjacent_atoms, incident_bonds, bond.end, bond.start, bond.id);
        }

        const component_of = allocator.alloc(MoleculeId, atoms.len) catch return error.OutOfMemory;
        errdefer allocator.free(component_of);
        @memset(component_of, .invalid);
        const component_offsets = allocator.alloc(u32, atoms.len + 1) catch return error.OutOfMemory;
        errdefer allocator.free(component_offsets);
        const component_atoms = allocator.alloc(AtomId, visibleAtomCount(atoms)) catch return error.OutOfMemory;
        errdefer allocator.free(component_atoms);

        var component_count: u32 = 0;
        var queued: u32 = 0;
        for (atoms, 0..) |atom, atom_index| {
            if (atom.hidden or component_of[atom_index] != .invalid) continue;
            component_offsets[component_count] = queued;
            const component_id = MoleculeId.fromIndex(component_count);
            component_count += 1;
            component_atoms[queued] = AtomId.fromIndex(@intCast(atom_index));
            component_of[atom_index] = component_id;
            queued += 1;

            var head = component_offsets[component_id.index()];
            while (head < queued) : (head += 1) {
                const current = component_atoms[head];
                const start = adjacency_offsets[current.index()];
                const end = adjacency_offsets[current.index() + 1];
                for (adjacent_atoms[start..end]) |neighbor| {
                    if (component_of[neighbor.index()] != .invalid) continue;
                    component_of[neighbor.index()] = component_id;
                    component_atoms[queued] = neighbor;
                    queued += 1;
                }
            }
            std.mem.sort(
                AtomId,
                component_atoms[component_offsets[component_id.index()]..queued],
                {},
                atomIdLessThan,
            );
        }
        component_offsets[component_count] = queued;

        return .{
            .allocator = allocator,
            .adjacency_offsets = adjacency_offsets,
            .adjacent_atoms = adjacent_atoms,
            .incident_bonds = incident_bonds,
            .structural_bonds = structural_bonds,
            .component_of = component_of,
            .component_offsets = component_offsets,
            .component_atoms = component_atoms,
            .component_count = component_count,
        };
    }

    pub fn deinit(self: *Graph) void {
        self.allocator.free(self.component_atoms);
        self.allocator.free(self.component_offsets);
        self.allocator.free(self.component_of);
        self.allocator.free(self.structural_bonds);
        self.allocator.free(self.incident_bonds);
        self.allocator.free(self.adjacent_atoms);
        self.allocator.free(self.adjacency_offsets);
        self.* = undefined;
    }

    pub fn neighbors(self: Graph, atom: AtomId) []const AtomId {
        const index = atom.index();
        std.debug.assert(index < self.adjacency_offsets.len - 1);
        return self.adjacent_atoms[self.adjacency_offsets[index]..self.adjacency_offsets[index + 1]];
    }

    pub fn incidentBonds(self: Graph, atom: AtomId) []const BondId {
        const index = atom.index();
        std.debug.assert(index < self.adjacency_offsets.len - 1);
        return self.incident_bonds[self.adjacency_offsets[index]..self.adjacency_offsets[index + 1]];
    }

    pub fn degree(self: Graph, atom: AtomId) u32 {
        const index = atom.index();
        std.debug.assert(index < self.adjacency_offsets.len - 1);
        return self.adjacency_offsets[index + 1] - self.adjacency_offsets[index];
    }

    pub fn structuralBonds(self: Graph) []const BondId {
        return self.structural_bonds;
    }

    pub fn component(self: Graph, atom: AtomId) ?MoleculeId {
        if (atom.index() >= self.component_of.len) return null;
        const result = self.component_of[atom.index()];
        return if (result == .invalid) null else result;
    }

    pub fn componentMembers(self: Graph, component_id: MoleculeId) []const AtomId {
        const index = component_id.index();
        std.debug.assert(component_id.isValid());
        std.debug.assert(index < self.component_count);
        return self.component_atoms[self.component_offsets[index]..self.component_offsets[index + 1]];
    }

    /// Breadth-first traversal while treating one atom as already visited,
    /// matching sketcherMinimizerAtom::getSubmolecule's cut-edge primitive.
    /// The caller owns the result and must free it with `allocator`.
    pub fn reachableExcluding(
        self: Graph,
        allocator: std.mem.Allocator,
        start: AtomId,
        excluded: AtomId,
    ) core.errors.Error![]AtomId {
        if (start.index() >= self.component_of.len or
            excluded.index() >= self.component_of.len or
            self.component(start) == null or
            self.component(excluded) == null)
        {
            return error.InvalidAtomIndex;
        }
        const visited = allocator.alloc(bool, self.component_of.len) catch return error.OutOfMemory;
        defer allocator.free(visited);
        @memset(visited, false);
        visited[excluded.index()] = true;
        visited[start.index()] = true;

        const result = allocator.alloc(AtomId, self.component_of.len) catch return error.OutOfMemory;
        errdefer allocator.free(result);
        var count: usize = 1;
        result[0] = start;
        var head: usize = 0;
        while (head < count) : (head += 1) {
            for (self.neighbors(result[head])) |neighbor| {
                if (visited[neighbor.index()]) continue;
                visited[neighbor.index()] = true;
                result[count] = neighbor;
                count += 1;
            }
        }
        if (count == result.len) return result;
        return allocator.realloc(result, count) catch return error.OutOfMemory;
    }
};

fn isStructuralBond(atoms: []const model.Atom, bond: model.Bond) bool {
    return !bond.skip and
        bond.effective_order != .zero and
        !atoms[bond.start.index()].hidden and
        !atoms[bond.end.index()].hidden;
}

fn appendEdge(
    cursors: []u32,
    adjacent_atoms: []AtomId,
    incident_bonds: []BondId,
    from: AtomId,
    to: AtomId,
    bond: BondId,
) void {
    const destination = cursors[from.index()];
    adjacent_atoms[destination] = to;
    incident_bonds[destination] = bond;
    cursors[from.index()] += 1;
}

fn visibleAtomCount(atoms: []const model.Atom) usize {
    var count: usize = 0;
    for (atoms) |atom| count += @intFromBool(!atom.hidden);
    return count;
}

fn atomIdLessThan(_: void, left: AtomId, right: AtomId) bool {
    return left.index() < right.index();
}

fn makeAtom(index: u32) model.Atom {
    return .{
        .id = AtomId.fromIndex(index),
        .input_index = index,
        .atomic_number = .carbon,
    };
}

fn makeBond(index: u32, start: u32, end: u32) model.Bond {
    return .{
        .id = BondId.fromIndex(index),
        .input_index = index,
        .start = AtomId.fromIndex(start),
        .end = AtomId.fromIndex(end),
        .input_order = .single,
        .effective_order = .single,
    };
}

test "adjacency filters nonstructural elements and retains paired bond order" {
    var atoms = [_]model.Atom{ makeAtom(0), makeAtom(1), makeAtom(2), makeAtom(3), makeAtom(4) };
    atoms[4].hidden = true;
    var input_bonds = [_]model.Bond{
        makeBond(0, 0, 1),
        makeBond(1, 0, 2),
        makeBond(2, 0, 3),
        makeBond(3, 1, 4),
        makeBond(4, 1, 2),
    };
    input_bonds[1].effective_order = .zero;
    input_bonds[2].skip = true;

    var graph = try Graph.init(std.testing.allocator, &atoms, &input_bonds);
    defer graph.deinit();

    try std.testing.expectEqualSlices(AtomId, &.{AtomId.fromIndex(1)}, graph.neighbors(AtomId.fromIndex(0)));
    try std.testing.expectEqualSlices(AtomId, &.{ AtomId.fromIndex(0), AtomId.fromIndex(2) }, graph.neighbors(AtomId.fromIndex(1)));
    try std.testing.expectEqualSlices(BondId, &.{ BondId.fromIndex(0), BondId.fromIndex(4) }, graph.incidentBonds(AtomId.fromIndex(1)));
    try std.testing.expectEqualSlices(BondId, &.{ BondId.fromIndex(0), BondId.fromIndex(4) }, graph.structuralBonds());
    try std.testing.expectEqual(@as(u32, 2), graph.component_count);
    try std.testing.expectEqualSlices(AtomId, &.{ AtomId.fromIndex(0), AtomId.fromIndex(1), AtomId.fromIndex(2) }, graph.componentMembers(MoleculeId.fromIndex(0)));
    try std.testing.expectEqualSlices(AtomId, &.{AtomId.fromIndex(3)}, graph.componentMembers(MoleculeId.fromIndex(1)));
    try std.testing.expect(graph.component(AtomId.fromIndex(4)) == null);
}

test "parallel bonds remain distinct and component traversal is deterministic" {
    const atoms = [_]model.Atom{ makeAtom(0), makeAtom(1), makeAtom(2) };
    const input_bonds = [_]model.Bond{ makeBond(0, 0, 1), makeBond(1, 0, 1), makeBond(2, 1, 2) };
    var graph = try Graph.init(std.testing.allocator, &atoms, &input_bonds);
    defer graph.deinit();

    try std.testing.expectEqual(@as(u32, 2), graph.degree(AtomId.fromIndex(0)));
    try std.testing.expectEqualSlices(BondId, &.{ BondId.fromIndex(0), BondId.fromIndex(1) }, graph.incidentBonds(AtomId.fromIndex(0)));
    try std.testing.expectEqualSlices(AtomId, &.{ AtomId.fromIndex(0), AtomId.fromIndex(1), AtomId.fromIndex(2) }, graph.componentMembers(MoleculeId.fromIndex(0)));

    const left = try graph.reachableExcluding(std.testing.allocator, AtomId.fromIndex(0), AtomId.fromIndex(1));
    defer std.testing.allocator.free(left);
    try std.testing.expectEqualSlices(AtomId, &.{AtomId.fromIndex(0)}, left);
    const right = try graph.reachableExcluding(std.testing.allocator, AtomId.fromIndex(2), AtomId.fromIndex(1));
    defer std.testing.allocator.free(right);
    try std.testing.expectEqualSlices(AtomId, &.{AtomId.fromIndex(2)}, right);
}

test "cyclic and disconnected components repeat in stable traversal order" {
    const atoms = [_]model.Atom{ makeAtom(0), makeAtom(1), makeAtom(2), makeAtom(3), makeAtom(4), makeAtom(5) };
    const input_bonds = [_]model.Bond{
        makeBond(0, 2, 0),
        makeBond(1, 4, 2),
        makeBond(2, 0, 4),
        makeBond(3, 1, 3),
        makeBond(4, 3, 5),
    };
    var first = try Graph.init(std.testing.allocator, &atoms, &input_bonds);
    defer first.deinit();
    var second = try Graph.init(std.testing.allocator, &atoms, &input_bonds);
    defer second.deinit();

    try std.testing.expectEqual(@as(u32, 2), first.component_count);
    try std.testing.expectEqualSlices(AtomId, &.{ AtomId.fromIndex(0), AtomId.fromIndex(2), AtomId.fromIndex(4) }, first.componentMembers(MoleculeId.fromIndex(0)));
    try std.testing.expectEqualSlices(AtomId, &.{ AtomId.fromIndex(1), AtomId.fromIndex(3), AtomId.fromIndex(5) }, first.componentMembers(MoleculeId.fromIndex(1)));
    try std.testing.expectEqualSlices(u32, first.adjacency_offsets, second.adjacency_offsets);
    try std.testing.expectEqualSlices(AtomId, first.adjacent_atoms, second.adjacent_atoms);
    try std.testing.expectEqualSlices(BondId, first.incident_bonds, second.incident_bonds);
    try std.testing.expectEqualSlices(BondId, first.structural_bonds, second.structural_bonds);
    try std.testing.expectEqualSlices(MoleculeId, first.component_of, second.component_of);
    try std.testing.expectEqualSlices(AtomId, first.component_atoms, second.component_atoms);
}

test "adjacency and component invariants hold across path and cycle families" {
    for (1..65) |atom_count| {
        const bond_count = if (atom_count > 2) atom_count else atom_count - 1;
        const atoms = try std.testing.allocator.alloc(model.Atom, atom_count);
        defer std.testing.allocator.free(atoms);
        const bonds = try std.testing.allocator.alloc(model.Bond, bond_count);
        defer std.testing.allocator.free(bonds);
        for (atoms, 0..) |*atom, index| atom.* = makeAtom(@intCast(index));
        for (bonds[0 .. atom_count - 1], 0..) |*bond, index| {
            bond.* = makeBond(@intCast(index), @intCast(index), @intCast(index + 1));
        }
        if (atom_count > 2) {
            bonds[bonds.len - 1] = makeBond(
                @intCast(bonds.len - 1),
                @intCast(atom_count - 1),
                0,
            );
        }

        var graph = try Graph.init(std.testing.allocator, atoms, bonds);
        defer graph.deinit();
        try std.testing.expectEqual(@as(u32, 1), graph.component_count);
        try std.testing.expectEqual(atom_count, graph.componentMembers(MoleculeId.fromIndex(0)).len);

        var degree_sum: usize = 0;
        for (atoms, 0..) |_, atom_index| {
            const atom = AtomId.fromIndex(@intCast(atom_index));
            degree_sum += graph.degree(atom);
            for (graph.neighbors(atom), graph.incidentBonds(atom)) |neighbor, bond_id| {
                const bond = bonds[bond_id.index()];
                try std.testing.expect(
                    (bond.start == atom and bond.end == neighbor) or
                        (bond.end == atom and bond.start == neighbor),
                );
            }
        }
        try std.testing.expectEqual(bond_count * 2, degree_sum);
    }
}

test "empty and malformed internal graphs are handled without assertions" {
    var empty = try Graph.init(std.testing.allocator, &.{}, &.{});
    defer empty.deinit();
    try std.testing.expectEqual(@as(u32, 0), empty.component_count);

    const atoms = [_]model.Atom{makeAtom(0)};
    try std.testing.expectError(error.InvalidAtomIndex, Graph.init(std.testing.allocator, &atoms, &.{makeBond(0, 0, 0)}));
    try std.testing.expectError(error.InvalidAtomIndex, Graph.init(std.testing.allocator, &atoms, &.{makeBond(0, 0, 1)}));
    var wrong_atom = makeAtom(1);
    try std.testing.expectError(error.InvalidMapping, Graph.init(std.testing.allocator, @as(*const [1]model.Atom, &wrong_atom), &.{}));

    const two_atoms = [_]model.Atom{ makeAtom(0), makeAtom(1) };
    const wrong_bond = makeBond(1, 0, 1);
    try std.testing.expectError(
        error.InvalidMapping,
        Graph.init(std.testing.allocator, &two_atoms, @as(*const [1]model.Bond, &wrong_bond)),
    );
    var graph = try Graph.init(std.testing.allocator, &two_atoms, &.{});
    defer graph.deinit();
    try std.testing.expectError(
        error.InvalidAtomIndex,
        graph.reachableExcluding(std.testing.allocator, AtomId.fromIndex(9), AtomId.fromIndex(0)),
    );
    try std.testing.expectError(
        error.InvalidAtomIndex,
        graph.reachableExcluding(std.testing.allocator, AtomId.fromIndex(0), AtomId.fromIndex(9)),
    );
}

fn constructAndDiscard(allocator: std.mem.Allocator) !void {
    const atoms = [_]model.Atom{ makeAtom(0), makeAtom(1), makeAtom(2), makeAtom(3) };
    const input_bonds = [_]model.Bond{ makeBond(0, 0, 1), makeBond(1, 1, 2), makeBond(2, 2, 3) };
    var graph = try Graph.init(allocator, &atoms, &input_bonds);
    graph.deinit();
}

test "graph construction reports and cleans up every allocation failure" {
    var counting_allocator = std.testing.FailingAllocator.init(std.testing.allocator, .{});
    try constructAndDiscard(counting_allocator.allocator());
    // Seven retained graph arrays plus the temporary adjacency cursor array.
    try std.testing.expectEqual(@as(usize, 8), counting_allocator.alloc_index);
    try std.testing.checkAllAllocationFailures(std.testing.allocator, constructAndDiscard, .{});
}

fn traverseAndDiscard(allocator: std.mem.Allocator, graph: Graph) !void {
    const reached = try graph.reachableExcluding(allocator, AtomId.fromIndex(0), AtomId.fromIndex(1));
    allocator.free(reached);
}

test "cut traversal reports and cleans up every allocation failure" {
    const atoms = [_]model.Atom{ makeAtom(0), makeAtom(1), makeAtom(2) };
    const input_bonds = [_]model.Bond{ makeBond(0, 0, 1), makeBond(1, 1, 2) };
    var graph = try Graph.init(std.testing.allocator, &atoms, &input_bonds);
    defer graph.deinit();
    var counting_allocator = std.testing.FailingAllocator.init(std.testing.allocator, .{});
    try traverseAndDiscard(counting_allocator.allocator(), graph);
    // The visited bitmap and caller-owned result are independently fallible.
    try std.testing.expectEqual(@as(usize, 2), counting_allocator.alloc_index);
    try std.testing.checkAllAllocationFailures(std.testing.allocator, traverseAndDiscard, .{graph});
}
