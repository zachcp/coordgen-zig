const std = @import("std");
const core = @import("core");
const model = @import("model");

pub const canonical = @import("topology/canonical.zig");
pub const prepare = @import("topology/prepare.zig");
pub const rings = @import("topology/rings.zig");
pub const stereo = @import("topology/stereo.zig");

pub const AtomId = core.ids.AtomId;
pub const BondId = core.ids.BondId;
pub const MoleculeId = core.ids.MoleculeId;

test {
    _ = canonical;
    _ = prepare;
    _ = rings;
    _ = stereo;
}

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

/// Complete preparation seam consumed by the generator: canonical owned model
/// storage plus structural adjacency and connected components built from that
/// exact storage. Deinitialization is the reverse construction order.
pub const PreparedGraph = struct {
    working: model.WorkingGraph,
    graph: Graph,
    rings: RingMembership,
    /// Probe-compatible scores concatenated in component order. Upstream's
    /// morganScores returns no vector for singleton components, so neither do
    /// we; this is deliberately not one value per active atom.
    component_morgan_scores: []u32,

    pub fn deinit(self: *PreparedGraph) void {
        self.working.allocator.free(self.component_morgan_scores);
        self.rings.deinit();
        self.graph.deinit();
        self.working.deinit();
        self.* = undefined;
    }
};

pub fn prepareInput(allocator: std.mem.Allocator, input: anytype) core.errors.Error!PreparedGraph {
    var normalized = try prepare.init(allocator, input);
    errdefer normalized.deinit();
    var graph = try Graph.init(allocator, normalized.working.atoms, normalized.working.bonds);
    errdefer graph.deinit();
    var ring_membership = try RingMembership.init(allocator, graph, normalized.working.bonds);
    errdefer ring_membership.deinit();
    const component_morgan_scores = try componentMorganScores(allocator, graph, normalized.working.bonds);
    return .{
        .working = normalized.working,
        .graph = graph,
        .rings = ring_membership,
        .component_morgan_scores = component_morgan_scores,
    };
}

pub const RingMembership = struct {
    allocator: std.mem.Allocator,
    rings: []model.Ring,
    ring_atoms: []AtomId,
    ring_bonds: []BondId,
    atom_ring_offsets: []u32,
    atom_rings: []core.ids.RingId,
    bond_ring_offsets: []u32,
    bond_rings: []core.ids.RingId,

    pub fn init(
        allocator: std.mem.Allocator,
        graph: Graph,
        bonds: []const model.Bond,
    ) core.errors.Error!RingMembership {
        var ring_list: std.ArrayList(model.Ring) = .empty;
        defer ring_list.deinit(allocator);
        var ring_bond_list: std.ArrayList(BondId) = .empty;
        defer ring_bond_list.deinit(allocator);
        const visited = allocator.alloc(bool, bonds.len) catch return error.OutOfMemory;
        defer allocator.free(visited);
        const parent = allocator.alloc(BondId, bonds.len) catch return error.OutOfMemory;
        defer allocator.free(parent);
        const parent_at_start = allocator.alloc(bool, bonds.len) catch return error.OutOfMemory;
        defer allocator.free(parent_at_start);
        const queue = allocator.alloc(BondId, graph.structural_bonds.len) catch return error.OutOfMemory;
        defer allocator.free(queue);
        const candidate = allocator.alloc(BondId, graph.structural_bonds.len) catch return error.OutOfMemory;
        defer allocator.free(candidate);

        for (graph.structuralBonds()) |root| {
            @memset(visited, false);
            @memset(parent, .invalid);
            @memset(parent_at_start, true);
            visited[root.index()] = true;
            queue[0] = root;
            var head: usize = 0;
            var tail: usize = 1;
            var closed = false;
            while (head < tail and !closed) : (head += 1) {
                const last = queue[head];
                const last_bond = bonds[last.index()];
                const pivot = if (parent_at_start[last.index()]) last_bond.end else last_bond.start;
                for (graph.incidentBonds(pivot)) |next| {
                    if (next == last) continue;
                    if (visited[next.index()]) {
                        if (next == root) {
                            var candidate_count: usize = 0;
                            var cursor = last;
                            while (cursor != .invalid) {
                                candidate[candidate_count] = cursor;
                                candidate_count += 1;
                                cursor = parent[cursor.index()];
                            }
                            if (!containsRing(ring_list.items, ring_bond_list.items, candidate[0..candidate_count])) {
                                const ring_id = core.ids.RingId.fromIndex(@intCast(ring_list.items.len));
                                ring_list.append(allocator, .{
                                    .id = ring_id,
                                    .atom_start = 0,
                                    .atom_count = 0,
                                    .bond_start = @intCast(ring_bond_list.items.len),
                                    .bond_count = @intCast(candidate_count),
                                }) catch return error.OutOfMemory;
                                ring_bond_list.appendSlice(allocator, candidate[0..candidate_count]) catch return error.OutOfMemory;
                            }
                            closed = true;
                        }
                    } else {
                        if (bonds[next.index()].end == pivot) parent_at_start[next.index()] = false;
                        parent[next.index()] = last;
                        visited[next.index()] = true;
                        queue[tail] = next;
                        tail += 1;
                    }
                }
            }
        }

        var ring_atom_list: std.ArrayList(AtomId) = .empty;
        defer ring_atom_list.deinit(allocator);
        for (ring_list.items) |*ring| {
            ring.atom_start = @intCast(ring_atom_list.items.len);
            const member_bonds = ring_bond_list.items[ring.bond_start..][0..ring.bond_count];
            for (0..graph.component_of.len) |atom_index| {
                const atom = AtomId.fromIndex(@intCast(atom_index));
                if (graph.component(atom) == null or !ringContainsAtom(member_bonds, bonds, atom)) continue;
                ring_atom_list.append(allocator, atom) catch return error.OutOfMemory;
            }
            ring.atom_count = @intCast(ring_atom_list.items.len - ring.atom_start);
        }

        const ring_records = ring_list.toOwnedSlice(allocator) catch return error.OutOfMemory;
        errdefer allocator.free(ring_records);
        const ring_bonds = ring_bond_list.toOwnedSlice(allocator) catch return error.OutOfMemory;
        errdefer allocator.free(ring_bonds);
        const ring_atoms = ring_atom_list.toOwnedSlice(allocator) catch return error.OutOfMemory;
        errdefer allocator.free(ring_atoms);
        const atom_membership = try buildMembership(allocator, graph.component_of.len, ring_records, ring_atoms, true);
        errdefer {
            allocator.free(atom_membership.ids);
            allocator.free(atom_membership.offsets);
        }
        const bond_membership = try buildMembership(allocator, bonds.len, ring_records, ring_bonds, false);
        return .{
            .allocator = allocator,
            .rings = ring_records,
            .ring_atoms = ring_atoms,
            .ring_bonds = ring_bonds,
            .atom_ring_offsets = atom_membership.offsets,
            .atom_rings = atom_membership.ids,
            .bond_ring_offsets = bond_membership.offsets,
            .bond_rings = bond_membership.ids,
        };
    }

    pub fn deinit(self: *RingMembership) void {
        self.allocator.free(self.bond_rings);
        self.allocator.free(self.bond_ring_offsets);
        self.allocator.free(self.atom_rings);
        self.allocator.free(self.atom_ring_offsets);
        self.allocator.free(self.ring_bonds);
        self.allocator.free(self.ring_atoms);
        self.allocator.free(self.rings);
        self.* = undefined;
    }

    pub fn atoms(self: RingMembership, ring: core.ids.RingId) []const AtomId {
        const record = self.rings[ring.index()];
        return self.ring_atoms[record.atom_start..][0..record.atom_count];
    }

    pub fn ringBonds(self: RingMembership, ring: core.ids.RingId) []const BondId {
        const record = self.rings[ring.index()];
        return self.ring_bonds[record.bond_start..][0..record.bond_count];
    }

    pub fn atomRings(self: RingMembership, atom: AtomId) []const core.ids.RingId {
        return self.atom_rings[self.atom_ring_offsets[atom.index()]..self.atom_ring_offsets[atom.index() + 1]];
    }

    pub fn bondRings(self: RingMembership, bond: BondId) []const core.ids.RingId {
        return self.bond_rings[self.bond_ring_offsets[bond.index()]..self.bond_ring_offsets[bond.index() + 1]];
    }
};

const Membership = struct { offsets: []u32, ids: []core.ids.RingId };

fn buildMembership(
    allocator: std.mem.Allocator,
    item_count: usize,
    ring_records: []const model.Ring,
    members: anytype,
    comptime atoms: bool,
) core.errors.Error!Membership {
    const offsets = allocator.alloc(u32, item_count + 1) catch return error.OutOfMemory;
    errdefer allocator.free(offsets);
    @memset(offsets, 0);
    for (ring_records) |ring| {
        const start = if (atoms) ring.atom_start else ring.bond_start;
        const count = if (atoms) ring.atom_count else ring.bond_count;
        for (members[start..][0..count]) |member| offsets[member.index() + 1] += 1;
    }
    for (offsets[1..], offsets[0 .. offsets.len - 1]) |*value, previous| value.* += previous;
    const ids = allocator.alloc(core.ids.RingId, offsets[item_count]) catch return error.OutOfMemory;
    errdefer allocator.free(ids);
    const cursors = allocator.dupe(u32, offsets[0..item_count]) catch return error.OutOfMemory;
    defer allocator.free(cursors);
    for (ring_records) |ring| {
        const start = if (atoms) ring.atom_start else ring.bond_start;
        const count = if (atoms) ring.atom_count else ring.bond_count;
        for (members[start..][0..count]) |member| {
            ids[cursors[member.index()]] = ring.id;
            cursors[member.index()] += 1;
        }
    }
    return .{ .offsets = offsets, .ids = ids };
}

fn containsRing(ring_records: []const model.Ring, ring_bonds: []const BondId, candidate: []const BondId) bool {
    for (ring_records) |ring| {
        if (ring.bond_count != candidate.len) continue;
        const existing = ring_bonds[ring.bond_start..][0..ring.bond_count];
        var same = true;
        for (candidate) |bond| {
            if (std.mem.indexOfScalar(BondId, existing, bond) == null) {
                same = false;
                break;
            }
        }
        if (same) return true;
    }
    return false;
}

fn ringContainsAtom(ring_bonds: []const BondId, bonds: []const model.Bond, atom: AtomId) bool {
    for (ring_bonds) |bond_id| {
        const bond = bonds[bond_id.index()];
        if (bond.start == atom or bond.end == atom) return true;
    }
    return false;
}

fn componentMorganScores(
    allocator: std.mem.Allocator,
    graph: Graph,
    bonds: []const model.Bond,
) core.errors.Error![]u32 {
    var score_count: usize = 0;
    var largest_component: usize = 0;
    for (0..graph.component_count) |component_index| {
        const members = graph.componentMembers(MoleculeId.fromIndex(@intCast(component_index)));
        if (members.len >= 2) score_count += members.len;
        largest_component = @max(largest_component, members.len);
    }
    const scores = allocator.alloc(u32, score_count) catch return error.OutOfMemory;
    errdefer allocator.free(scores);
    if (score_count == 0) return scores;
    const cumulative = allocator.alloc(u32, largest_component) catch return error.OutOfMemory;
    defer allocator.free(cumulative);
    const ordered = allocator.alloc(u32, largest_component) catch return error.OutOfMemory;
    defer allocator.free(ordered);

    var score_offset: usize = 0;
    for (0..graph.component_count) |component_index| {
        const component_id = MoleculeId.fromIndex(@intCast(component_index));
        const members = graph.componentMembers(component_id);
        if (members.len < 2) continue;
        const current = scores[score_offset..][0..members.len];
        @memset(current, 1);
        @memset(cumulative[0..members.len], 0);
        var old_ties = members.len;
        while (true) {
            for (graph.structuralBonds()) |bond_id| {
                const bond = bonds[bond_id.index()];
                if (graph.component(bond.start).? != component_id) continue;
                const start = memberIndex(members, bond.start);
                const end = memberIndex(members, bond.end);
                cumulative[start] +%= current[end];
                cumulative[end] +%= current[start];
            }
            @memcpy(ordered[0..members.len], cumulative[0..members.len]);
            std.mem.sort(u32, ordered[0..members.len], {}, std.sort.asc(u32));
            var new_ties: usize = 0;
            for (ordered[1..members.len], ordered[0 .. members.len - 1]) |score, previous| {
                new_ties += @intFromBool(score == previous);
            }
            if (new_ties >= old_ties) break;
            old_ties = new_ties;
            @memcpy(current, cumulative[0..members.len]);
        }
        score_offset += members.len;
    }
    return scores;
}

fn memberIndex(members: []const AtomId, atom: AtomId) usize {
    for (members, 0..) |member, index| if (member == atom) return index;
    unreachable;
}

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

test "component Morgan scores match probe semantics including singleton omission" {
    const atoms = [_]model.Atom{ makeAtom(0), makeAtom(1), makeAtom(2), makeAtom(3) };
    const input_bonds = [_]model.Bond{ makeBond(0, 0, 1), makeBond(1, 0, 2) };
    var graph = try Graph.init(std.testing.allocator, &atoms, &input_bonds);
    defer graph.deinit();
    const scores = try componentMorganScores(std.testing.allocator, graph, &input_bonds);
    defer std.testing.allocator.free(scores);
    try std.testing.expectEqualSlices(u32, &.{ 2, 1, 1 }, scores);
}

fn scoreAndDiscard(allocator: std.mem.Allocator, graph: Graph, bonds: []const model.Bond) !void {
    const scores = try componentMorganScores(allocator, graph, bonds);
    allocator.free(scores);
}

test "component Morgan scoring reports and cleans every allocation failure" {
    const atoms = [_]model.Atom{ makeAtom(0), makeAtom(1), makeAtom(2), makeAtom(3) };
    const input_bonds = [_]model.Bond{ makeBond(0, 0, 1), makeBond(1, 0, 2) };
    var graph = try Graph.init(std.testing.allocator, &atoms, &input_bonds);
    defer graph.deinit();
    try std.testing.checkAllAllocationFailures(
        std.testing.allocator,
        scoreAndDiscard,
        .{ graph, &input_bonds },
    );
}

fn prepareAndDiscard(
    allocator: std.mem.Allocator,
    atoms: []const prepare.TestAtom,
    bonds: []const prepare.TestBond,
) !void {
    var prepared = try prepareInput(allocator, prepare.TestInput{ .atoms = atoms, .bonds = bonds });
    prepared.deinit();
}

test "complete preparation seam matches pinned oracle maps, components, and Morgan ranks" {
    const atoms = [_]prepare.TestAtom{
        .{ .atomic_number = .iron },
        .{ .atomic_number = .carbon },
        .{ .atomic_number = .oxygen },
        .{ .atomic_number = .phosphorus },
        .{ .atomic_number = .carbon, .hidden = true },
    };
    const bonds = [_]prepare.TestBond{
        .{ .start = 0, .end = 1 },
        .{ .start = 0, .end = 2 },
        .{ .start = 1, .end = 3 },
        .{ .start = 3, .end = 4 },
    };
    var prepared = try prepareInput(std.testing.allocator, prepare.TestInput{ .atoms = &atoms, .bonds = &bonds });
    defer prepared.deinit();

    try std.testing.expectEqualSlices(u32, &.{ 0, 2, 3, 1 }, prepared.working.order.activeInternalToInput(4));
    try std.testing.expectEqual(@as(u32, 2), prepared.graph.component_count);
    try std.testing.expectEqualSlices(AtomId, &.{ AtomId.fromIndex(0), AtomId.fromIndex(1) }, prepared.graph.componentMembers(MoleculeId.fromIndex(0)));
    try std.testing.expectEqualSlices(AtomId, &.{ AtomId.fromIndex(2), AtomId.fromIndex(3) }, prepared.graph.componentMembers(MoleculeId.fromIndex(1)));
    try std.testing.expectEqualSlices(u32, &.{ 1, 1, 1, 1 }, prepared.component_morgan_scores);

    try std.testing.checkAllAllocationFailures(
        std.testing.allocator,
        prepareAndDiscard,
        .{ &atoms, &bonds },
    );
}

test "ring perception matches pinned oracle cycle and fused probes" {
    const cycle_atoms = [_]prepare.TestAtom{ .{}, .{ .atomic_number = .nitrogen }, .{ .atomic_number = .oxygen }, .{} };
    const cycle_bonds = [_]prepare.TestBond{
        .{ .start = 0, .end = 1 },
        .{ .start = 1, .end = 2 },
        .{ .start = 2, .end = 3 },
        .{ .start = 3, .end = 0 },
    };
    var cycle = try prepareInput(std.testing.allocator, prepare.TestInput{ .atoms = &cycle_atoms, .bonds = &cycle_bonds });
    defer cycle.deinit();
    try std.testing.expectEqual(@as(usize, 1), cycle.rings.rings.len);
    try std.testing.expectEqualSlices(AtomId, &.{ AtomId.fromIndex(0), AtomId.fromIndex(1), AtomId.fromIndex(2), AtomId.fromIndex(3) }, cycle.rings.atoms(core.ids.RingId.fromIndex(0)));
    try std.testing.expectEqualSlices(core.ids.RingId, &.{core.ids.RingId.fromIndex(0)}, cycle.rings.atomRings(AtomId.fromIndex(0)));

    const fused_atoms = [_]prepare.TestAtom{ .{}, .{}, .{}, .{}, .{}, .{} };
    const fused_bonds = [_]prepare.TestBond{
        .{ .start = 0, .end = 1 },
        .{ .start = 1, .end = 2 },
        .{ .start = 2, .end = 3 },
        .{ .start = 3, .end = 0 },
        .{ .start = 2, .end = 5 },
        .{ .start = 5, .end = 4 },
        .{ .start = 4, .end = 3 },
    };
    var fused = try prepareInput(std.testing.allocator, prepare.TestInput{ .atoms = &fused_atoms, .bonds = &fused_bonds });
    defer fused.deinit();
    try std.testing.expectEqualSlices(u32, &.{ 2, 3, 1, 5, 0, 4 }, fused.working.order.activeInternalToInput(6));
    try std.testing.expectEqual(@as(usize, 2), fused.rings.rings.len);
    try std.testing.expectEqualSlices(AtomId, &.{ AtomId.fromIndex(0), AtomId.fromIndex(1), AtomId.fromIndex(2), AtomId.fromIndex(4) }, fused.rings.atoms(core.ids.RingId.fromIndex(0)));
    try std.testing.expectEqualSlices(AtomId, &.{ AtomId.fromIndex(0), AtomId.fromIndex(1), AtomId.fromIndex(3), AtomId.fromIndex(5) }, fused.rings.atoms(core.ids.RingId.fromIndex(1)));
    try std.testing.expectEqual(@as(usize, 2), fused.rings.atomRings(AtomId.fromIndex(0)).len);
    var fusion = try rings.Analysis.init(std.testing.allocator, fused.rings, fused.working.atoms, fused.working.bonds);
    defer fusion.deinit();
    try std.testing.expectEqual(@as(usize, 1), fusion.fusedWith(core.ids.RingId.fromIndex(0)).len);
    try std.testing.expectEqualSlices(
        AtomId,
        &.{ AtomId.fromIndex(0), AtomId.fromIndex(1) },
        fusion.fusionAtoms(fusion.fusedWith(core.ids.RingId.fromIndex(0))[0]),
    );
    try std.testing.expect(!fusion.flags[0].aromatic);

    try std.testing.checkAllAllocationFailures(
        std.testing.allocator,
        prepareAndDiscard,
        .{ &fused_atoms, &fused_bonds },
    );
}

test "benzene ring chemistry matches pinned upstream heuristics" {
    const atoms = [_]prepare.TestAtom{ .{}, .{}, .{}, .{}, .{}, .{} };
    const bonds = [_]prepare.TestBond{
        .{ .start = 0, .end = 1, .order = .double },
        .{ .start = 1, .end = 2 },
        .{ .start = 2, .end = 3, .order = .double },
        .{ .start = 3, .end = 4 },
        .{ .start = 4, .end = 5, .order = .double },
        .{ .start = 5, .end = 0 },
    };
    var prepared = try prepareInput(std.testing.allocator, prepare.TestInput{ .atoms = &atoms, .bonds = &bonds });
    defer prepared.deinit();
    var analysis = try rings.Analysis.init(std.testing.allocator, prepared.rings, prepared.working.atoms, prepared.working.bonds);
    defer analysis.deinit();
    try std.testing.expectEqual(@as(usize, 1), analysis.flags.len);
    try std.testing.expect(analysis.flags[0].benzene);
    try std.testing.expect(analysis.flags[0].aromatic);
    try std.testing.expect(!analysis.flags[0].macrocycle);
}

test "CIP branch order and relative double-bond stereo produce absolute Z/E" {
    const atoms = [_]prepare.TestAtom{
        .{},
        .{},
        .{ .atomic_number = .oxygen },
        .{},
        .{ .atomic_number = .nitrogen },
        .{},
    };
    const bonds = [_]prepare.TestBond{
        .{ .start = 0, .end = 1, .order = .double, .stereo = .{ .value = .cis, .atom_a = 2, .atom_b = 4 } },
        .{ .start = 0, .end = 2 },
        .{ .start = 0, .end = 3 },
        .{ .start = 1, .end = 4 },
        .{ .start = 1, .end = 5 },
    };
    var prepared = try prepareInput(std.testing.allocator, prepare.TestInput{ .atoms = &atoms, .bonds = &bonds });
    defer prepared.deinit();
    var double_bond: BondId = .invalid;
    for (prepared.working.bonds) |bond| if (bond.input_index == 0) {
        double_bond = bond.id;
        break;
    };
    try std.testing.expect(double_bond.isValid());
    const bond = prepared.working.bonds[double_bond.index()];
    const start_neighbor = (try stereo.firstNeighbor(
        std.testing.allocator,
        prepared.working.atoms,
        prepared.working.bonds,
        prepared.graph,
        bond.start,
        bond.end,
    )).?;
    const end_neighbor = (try stereo.firstNeighbor(
        std.testing.allocator,
        prepared.working.atoms,
        prepared.working.bonds,
        prepared.graph,
        bond.end,
        bond.start,
    )).?;
    try std.testing.expectEqual(prepared.working.order.input_to_internal[2], start_neighbor.index());
    try std.testing.expectEqual(prepared.working.order.input_to_internal[4], end_neighbor.index());
    try std.testing.expectEqual(.z, try stereo.absoluteBondStereo(
        std.testing.allocator,
        prepared.working.atoms,
        prepared.working.bonds,
        prepared.graph,
        prepared.rings,
        double_bond,
    ));
}

fn tetrahedralFixture(stereo_value: core.chemistry.AtomStereo) struct {
    [5]model.Atom,
    [4]model.Bond,
} {
    var atoms = [_]model.Atom{ makeAtom(0), makeAtom(1), makeAtom(2), makeAtom(3), makeAtom(4) };
    atoms[0].stereo = stereo_value;
    atoms[0].stereo_atom_a = AtomId.fromIndex(4);
    atoms[0].stereo_atom_b = AtomId.fromIndex(3);
    atoms[0].stereo_looking_from = AtomId.fromIndex(1);
    atoms[1].atomic_number = .fluorine;
    atoms[2].atomic_number = .chlorine;
    atoms[3].atomic_number = .bromine;
    atoms[4].atomic_number = .iodine;
    atoms[0].coordinates = .{};
    atoms[1].coordinates = .{ .x = 1 };
    atoms[2].coordinates = .{ .y = 1 };
    atoms[3].coordinates = .{ .x = -1 };
    atoms[4].coordinates = .{ .y = -1 };
    const bonds = [_]model.Bond{
        makeBond(0, 0, 1),
        makeBond(1, 0, 2),
        makeBond(2, 0, 3),
        makeBond(3, 0, 4),
    };
    return .{ atoms, bonds };
}

test "relative atom stereo resolves to absolute CIP and mirrored input inverts" {
    var clockwise = tetrahedralFixture(.counter_clockwise);
    var graph = try Graph.init(std.testing.allocator, &clockwise[0], &clockwise[1]);
    defer graph.deinit();
    try std.testing.expectEqual(.r, try stereo.absoluteAtomStereo(
        std.testing.allocator,
        &clockwise[0],
        &clockwise[1],
        graph,
        AtomId.fromIndex(0),
    ));

    var mirrored = tetrahedralFixture(.clockwise);
    try std.testing.expectEqual(.s, try stereo.absoluteAtomStereo(
        std.testing.allocator,
        &mirrored[0],
        &mirrored[1],
        graph,
        AtomId.fromIndex(0),
    ));
}

test "absolute atom stereo writes deterministic displays and pseudo depth" {
    var fixture = tetrahedralFixture(.counter_clockwise);
    var graph = try Graph.init(std.testing.allocator, &fixture[0], &fixture[1]);
    defer graph.deinit();
    var membership = try RingMembership.init(std.testing.allocator, graph, &fixture[1]);
    defer membership.deinit();
    try stereo.writeAtomBondDisplays(std.testing.allocator, &fixture[0], &fixture[1], graph, membership);
    try std.testing.expectEqual(.r, fixture[0][0].stereo);
    var displayed: usize = 0;
    for (fixture[1]) |bond| displayed += @intFromBool(bond.display != .none);
    try std.testing.expectEqual(@as(usize, 2), displayed);

    const pseudo_z = try stereo.assignPseudoZ(std.testing.allocator, &fixture[0], &fixture[1], graph);
    defer std.testing.allocator.free(pseudo_z);
    try std.testing.expectEqual(@as(f32, 0), pseudo_z[0]);
    try std.testing.expectEqual(@as(f32, 1), @abs(pseudo_z[2]) + @abs(pseudo_z[1]) + @abs(pseudo_z[3]) + @abs(pseudo_z[4]) - 1);
}

fn resolveTetrahedralAndDiscard(
    allocator: std.mem.Allocator,
    atoms: []const model.Atom,
    bonds: []const model.Bond,
    graph: Graph,
) !void {
    _ = try stereo.absoluteAtomStereo(allocator, atoms, bonds, graph, AtomId.fromIndex(0));
}

fn writeTetrahedralAndDiscard(
    allocator: std.mem.Allocator,
    source_atoms: []const model.Atom,
    source_bonds: []const model.Bond,
    graph: Graph,
    membership: RingMembership,
) !void {
    const atoms = try allocator.dupe(model.Atom, source_atoms);
    defer allocator.free(atoms);
    const bonds = try allocator.dupe(model.Bond, source_bonds);
    defer allocator.free(bonds);
    try stereo.writeAtomBondDisplays(allocator, atoms, bonds, graph, membership);
}

test "malformed and allocation-failing atom stereo is rejected cleanly" {
    var fixture = tetrahedralFixture(.counter_clockwise);
    var graph = try Graph.init(std.testing.allocator, &fixture[0], &fixture[1]);
    defer graph.deinit();
    fixture[0][0].stereo_atom_b = fixture[0][0].stereo_atom_a;
    try std.testing.expectError(error.InvalidStereo, stereo.absoluteAtomStereo(
        std.testing.allocator,
        &fixture[0],
        &fixture[1],
        graph,
        AtomId.fromIndex(0),
    ));
    fixture[0][0].stereo_atom_b = AtomId.fromIndex(3);
    try std.testing.checkAllAllocationFailures(
        std.testing.allocator,
        resolveTetrahedralAndDiscard,
        .{ &fixture[0], &fixture[1], graph },
    );
    var membership = try RingMembership.init(std.testing.allocator, graph, &fixture[1]);
    defer membership.deinit();
    try std.testing.checkAllAllocationFailures(
        std.testing.allocator,
        writeTetrahedralAndDiscard,
        .{ &fixture[0], &fixture[1], graph, membership },
    );
}

test "peptide mirror scoring follows the pinned N-to-carbonyl preference" {
    var atoms = [_]model.Atom{ makeAtom(0), makeAtom(1), makeAtom(2), makeAtom(3) };
    atoms[1].atomic_number = .nitrogen;
    atoms[3].atomic_number = .oxygen;
    atoms[1].coordinates.x = 1;
    atoms[2].coordinates.x = -1;
    var bonds = [_]model.Bond{
        makeBond(0, 0, 1),
        makeBond(1, 0, 2),
        makeBond(2, 2, 3),
    };
    bonds[2].input_order = .double;
    bonds[2].effective_order = .double;
    var graph = try Graph.init(std.testing.allocator, &atoms, &bonds);
    defer graph.deinit();
    try std.testing.expectEqual(@as(f32, -100), stereo.peptideFlipScore(&atoms, &bonds, graph));
    atoms[1].coordinates.x = -1;
    atoms[2].coordinates.x = 1;
    try std.testing.expectEqual(@as(f32, 100), stereo.peptideFlipScore(&atoms, &bonds, graph));
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
