const std = @import("std");
const core = @import("core");
const model = @import("model");
const topology = @import("topology");

pub const Flags = packed struct {
    fixed: bool = false,
    constrained: bool = false,
    chain: bool = false,
    constrained_flip: bool = false,
};

pub const Fragment = struct {
    id: core.ids.FragmentId,
    component: core.ids.MoleculeId,
    parent: core.ids.FragmentId = .invalid,
    bond_to_parent: core.ids.BondId = .invalid,
    /// The two ends of `bond_to_parent`, resolved once by fragment membership
    /// so no caller has to work them out from the bond's stored direction.
    ///
    /// Upstream can read `_bondToParent->startAtom` directly because
    /// CoordgenFragmenter mutates the molecule to guarantee it is the parent's
    /// end (CoordgenFragmenter.cpp:406-422, with an assert). Native keeps
    /// bonds in canonical input order, which on the drug_like corpus is
    /// child-first, so that guarantee does not hold and reading `.start`
    /// selects the wrong atom - silently (cgz-jg4).
    attachment_atom: core.ids.AtomId = .invalid,
    parent_atom: core.ids.AtomId = .invalid,
    atom_start: u32,
    atom_count: u32,
    ring_count: u32,
    inter_bond_count: u32,
    flags: Flags,
};

pub const Fragmentation = struct {
    allocator: std.mem.Allocator,
    fragments: []Fragment,
    atoms: []core.ids.AtomId,
    atom_fragment: []core.ids.FragmentId,
    main_fragments: []core.ids.FragmentId,

    pub fn init(
        allocator: std.mem.Allocator,
        atoms: []const model.Atom,
        bonds: []const model.Bond,
        graph: anytype,
        membership: anytype,
    ) core.errors.Error!Fragmentation {
        const parents = allocator.alloc(u32, atoms.len) catch return error.OutOfMemory;
        defer allocator.free(parents);
        for (parents, 0..) |*parent, index| parent.* = @intCast(index);
        for (graph.structuralBonds()) |bond_id| {
            const bond = bonds[bond_id.index()];
            if (!isInterFragmentBond(graph, membership, bond)) join(parents, bond.start.index(), bond.end.index());
        }

        const atom_fragment = allocator.alloc(core.ids.FragmentId, atoms.len) catch return error.OutOfMemory;
        errdefer allocator.free(atom_fragment);
        @memset(atom_fragment, .invalid);
        const root_fragment = allocator.alloc(core.ids.FragmentId, atoms.len) catch return error.OutOfMemory;
        defer allocator.free(root_fragment);
        @memset(root_fragment, .invalid);
        var fragment_count: u32 = 0;
        // Fragment objects are born while upstream walks bonds, not while it
        // walks atoms. Preserve that first-appearance order after the union of
        // rigid groups has been established.
        for (graph.structuralBonds()) |bond_id| {
            const bond = bonds[bond_id.index()];
            const roots = [_]usize{ find(parents, bond.start.index()), find(parents, bond.end.index()) };
            for (roots) |root| {
                if (root_fragment[root].isValid()) continue;
                root_fragment[root] = core.ids.FragmentId.fromIndex(fragment_count);
                fragment_count += 1;
            }
        }
        for (atoms) |atom| {
            if (atom.hidden or graph.component(atom.id) == null) continue;
            const root = find(parents, atom.id.index());
            if (!root_fragment[root].isValid()) {
                root_fragment[root] = core.ids.FragmentId.fromIndex(fragment_count);
                fragment_count += 1;
            }
            atom_fragment[atom.id.index()] = root_fragment[root];
        }

        const records = allocator.alloc(Fragment, fragment_count) catch return error.OutOfMemory;
        errdefer allocator.free(records);
        const counts = allocator.alloc(u32, fragment_count) catch return error.OutOfMemory;
        defer allocator.free(counts);
        @memset(counts, 0);
        for (atom_fragment) |fragment| {
            if (fragment.isValid()) counts[fragment.index()] += 1;
        }
        var offset: u32 = 0;
        for (records, 0..) |*record, index| {
            record.* = .{
                .id = core.ids.FragmentId.fromIndex(@intCast(index)),
                .component = .invalid,
                .atom_start = offset,
                .atom_count = counts[index],
                .ring_count = 0,
                .inter_bond_count = 0,
                .flags = .{},
            };
            offset += counts[index];
            counts[index] = record.atom_start;
        }
        const fragment_atoms = allocator.alloc(core.ids.AtomId, offset) catch return error.OutOfMemory;
        errdefer allocator.free(fragment_atoms);
        for (atoms) |atom| {
            const fragment = atom_fragment[atom.id.index()];
            if (!fragment.isValid()) continue;
            fragment_atoms[counts[fragment.index()]] = atom.id;
            counts[fragment.index()] += 1;
        }

        for (records) |*record| {
            const fragment_members = fragment_atoms[record.atom_start..][0..record.atom_count];
            record.component = graph.component(fragment_members[0]).?;
            var constrained_count: u32 = 0;
            for (fragment_members) |atom_id| {
                const atom = atoms[atom_id.index()];
                record.flags.fixed = record.flags.fixed or atom.fixed;
                record.flags.constrained = record.flags.constrained or atom.constrained;
                constrained_count += @intFromBool(atom.constrained);
                record.ring_count += @intCast(membership.atomRings(atom_id).len);
            }
            // Every ring was counted once per member above.
            for (membership.rings) |ring| {
                if (atom_fragment[membership.atoms(ring.id)[0].index()] == record.id) {
                    record.ring_count -= ring.atom_count - 1;
                }
            }
            record.flags.constrained_flip = if (record.atom_count == 1) false else constrained_count > 1;
        }
        for (graph.structuralBonds()) |bond_id| {
            const bond = bonds[bond_id.index()];
            const start_fragment = atom_fragment[bond.start.index()];
            const end_fragment = atom_fragment[bond.end.index()];
            if (start_fragment != end_fragment) {
                records[start_fragment.index()].inter_bond_count += 1;
                records[end_fragment.index()].inter_bond_count += 1;
            }
        }
        for (records) |*record| record.flags.chain = isChain(record.*, fragment_atoms, atoms, bonds, graph, membership);

        const main_fragments = allocator.alloc(core.ids.FragmentId, graph.component_count) catch return error.OutOfMemory;
        errdefer allocator.free(main_fragments);
        @memset(main_fragments, .invalid);
        for (records) |record| {
            const slot = &main_fragments[record.component.index()];
            if (!slot.isValid() or hasPriority(record, records[slot.index()], fragment_atoms, atoms, bonds, graph)) slot.* = record.id;
        }
        try considerChains(allocator, records, atom_fragment, main_fragments, bonds, graph);
        try assignParents(allocator, records, atom_fragment, main_fragments, bonds, graph);
        for (records) |*record| {
            if (record.atom_count != 1 or !record.flags.constrained) continue;
            for (records) |child| {
                if (child.parent == record.id and child.flags.constrained) {
                    record.flags.constrained_flip = true;
                    break;
                }
            }
        }
        return .{
            .allocator = allocator,
            .fragments = records,
            .atoms = fragment_atoms,
            .atom_fragment = atom_fragment,
            .main_fragments = main_fragments,
        };
    }

    pub fn deinit(self: *Fragmentation) void {
        self.allocator.free(self.main_fragments);
        self.allocator.free(self.atom_fragment);
        self.allocator.free(self.atoms);
        self.allocator.free(self.fragments);
        self.* = undefined;
    }

    pub fn members(self: Fragmentation, fragment: core.ids.FragmentId) []const core.ids.AtomId {
        const record = self.fragments[fragment.index()];
        return self.atoms[record.atom_start..][0..record.atom_count];
    }
};

pub fn isInterFragmentBond(graph: anytype, membership: anytype, bond: model.Bond) bool {
    if (graph.degree(bond.start) == 1 or graph.degree(bond.end) == 1) return false;
    if (membership.bondRings(bond.id).len != 0) return false;
    if (bond.effective_order == .double) return false;
    return true;
}

fn find(parents: []u32, start: usize) usize {
    var cursor = start;
    while (parents[cursor] != cursor) cursor = parents[cursor];
    return cursor;
}

fn join(parents: []u32, left: usize, right: usize) void {
    const left_root = find(parents, left);
    const right_root = find(parents, right);
    if (left_root != right_root) parents[right_root] = @intCast(left_root);
}

fn isChain(record: Fragment, members: []const core.ids.AtomId, atoms: []const model.Atom, bonds: []const model.Bond, graph: anytype, membership: anytype) bool {
    if (record.atom_count > 3) return false;
    for (members[record.atom_start..][0..record.atom_count]) |atom| {
        if (graph.degree(atom) > 3 or membership.atomRings(atom).len != 0) return false;
        _ = atoms;
    }
    for (graph.structuralBonds()) |bond_id| {
        const bond = bonds[bond_id.index()];
        if (std.mem.indexOfScalar(core.ids.AtomId, members[record.atom_start..][0..record.atom_count], bond.start) != null and
            std.mem.indexOfScalar(core.ids.AtomId, members[record.atom_start..][0..record.atom_count], bond.end) != null and
            bond.effective_order == .triple) return false;
    }
    return true;
}

fn hasPriority(left: Fragment, right: Fragment, members: []const core.ids.AtomId, atoms: []const model.Atom, bonds: []const model.Bond, graph: anytype) bool {
    const left_values = priorityValues(left, members, atoms, bonds, graph);
    const right_values = priorityValues(right, members, atoms, bonds, graph);
    for (left_values, right_values) |left_value, right_value| {
        if (left_value != right_value) return left_value > right_value;
    }
    return false;
}

fn priorityValues(record: Fragment, members: []const core.ids.AtomId, atoms: []const model.Atom, bonds: []const model.Bond, graph: anytype) [8]u64 {
    var fixed: u64 = 0;
    var constrained: u64 = 0;
    var heavy: u64 = 0;
    var weight: u64 = 0;
    var doubles: u64 = 0;
    const fragment_atoms = members[record.atom_start..][0..record.atom_count];
    for (fragment_atoms) |atom_id| {
        const atom = atoms[atom_id.index()];
        fixed += @intFromBool(atom.fixed);
        constrained += @intFromBool(atom.constrained);
        // Upstream's historical countHeavyAtoms name is misleading: it
        // counts non-carbon atoms, including explicit hydrogens.
        heavy += @intFromBool(atom.atomic_number != .carbon);
        weight += @backingInt(atom.atomic_number);
        for (graph.incidentBonds(atom_id)) |bond_id| {
            if (bonds[bond_id.index()].effective_order == .double) doubles += 1;
        }
    }
    return .{ fixed, constrained, record.ring_count, record.atom_count, record.inter_bond_count, heavy, weight, doubles / 2 };
}

fn considerChains(allocator: std.mem.Allocator, records: []const Fragment, atom_fragment: []const core.ids.FragmentId, mains: []core.ids.FragmentId, bonds: []const model.Bond, graph: anytype) core.errors.Error!void {
    const queue = allocator.alloc(core.ids.FragmentId, records.len) catch return error.OutOfMemory;
    defer allocator.free(queue);
    const visited = allocator.alloc(bool, records.len) catch return error.OutOfMemory;
    defer allocator.free(visited);
    for (mains, 0..) |*main, component_index| {
        if (!main.isValid()) continue;
        var constrained = false;
        for (records) |record| if (record.component.index() == component_index) {
            constrained = constrained or record.flags.fixed or record.flags.constrained;
        };
        if (constrained) continue;
        var longest: usize = 0;
        var longest_start = main.*;
        for (records) |record| {
            if (record.component.index() != component_index or !record.flags.chain or
                chainNeighborCount(record.id, records, atom_fragment, bonds, graph) > 1) continue;
            @memset(visited, false);
            visited[record.id.index()] = true;
            queue[0] = record.id;
            var head: usize = 0;
            var tail: usize = 1;
            while (head < tail) : (head += 1) {
                const current = queue[head];
                for (graph.structuralBonds()) |bond_id| {
                    const bond = bonds[bond_id.index()];
                    const start = atom_fragment[bond.start.index()];
                    const end = atom_fragment[bond.end.index()];
                    const next = if (start == current) end else if (end == current) start else continue;
                    if (next == current or visited[next.index()] or !records[next.index()].flags.chain) continue;
                    visited[next.index()] = true;
                    queue[tail] = next;
                    tail += 1;
                }
            }
            if (tail > longest) {
                longest = tail;
                longest_start = record.id;
            }
        }
        if (longest >= acceptableChainLength(records[main.index()].ring_count)) main.* = longest_start;
    }
}

fn chainNeighborCount(fragment: core.ids.FragmentId, records: []const Fragment, atom_fragment: []const core.ids.FragmentId, bonds: []const model.Bond, graph: anytype) usize {
    var count: usize = 0;
    for (graph.structuralBonds()) |bond_id| {
        const bond = bonds[bond_id.index()];
        const start = atom_fragment[bond.start.index()];
        const end = atom_fragment[bond.end.index()];
        const other = if (start == fragment) end else if (end == fragment) start else continue;
        if (other != fragment and records[other.index()].flags.chain) count += 1;
    }
    return count;
}

fn acceptableChainLength(ring_count: u32) usize {
    return switch (ring_count) {
        0 => 1,
        1 => 5,
        2 => 8,
        3 => 10,
        else => 12,
    };
}

fn assignParents(allocator: std.mem.Allocator, records: []Fragment, atom_fragment: []const core.ids.FragmentId, mains: []const core.ids.FragmentId, bonds: []const model.Bond, graph: anytype) core.errors.Error!void {
    const queue = allocator.alloc(core.ids.FragmentId, records.len) catch return error.OutOfMemory;
    defer allocator.free(queue);
    for (mains) |main| {
        if (!main.isValid()) continue;
        var head: usize = 0;
        var tail: usize = 1;
        queue[0] = main;
        while (head < tail) : (head += 1) {
            const current = queue[head];
            for (graph.structuralBonds()) |bond_id| {
                const bond = bonds[bond_id.index()];
                const start = atom_fragment[bond.start.index()];
                const end = atom_fragment[bond.end.index()];
                if (start == end) continue;
                const child = if (start == current) end else if (end == current) start else continue;
                if (child == main or records[child.index()].parent.isValid()) continue;
                records[child.index()].parent = current;
                records[child.index()].bond_to_parent = bond_id;
                // Orientation is known here - `child` is the fragment being
                // attached - so it is resolved once rather than rediscovered
                // at every use site.
                const child_is_start = start == child;
                records[child.index()].attachment_atom = if (child_is_start) bond.start else bond.end;
                records[child.index()].parent_atom = if (child_is_start) bond.end else bond.start;
                queue[tail] = child;
                tail += 1;
            }
        }
    }
}

fn testAtom(index: u32) model.Atom {
    return .{
        .id = core.ids.AtomId.fromIndex(index),
        .input_index = index,
        .atomic_number = .carbon,
    };
}

fn testBond(index: u32, start: u32, end: u32) model.Bond {
    return .{
        .id = core.ids.BondId.fromIndex(index),
        .input_index = index,
        .start = core.ids.AtomId.fromIndex(start),
        .end = core.ids.AtomId.fromIndex(end),
        .input_order = .single,
        .effective_order = .single,
    };
}

fn fragmentAndDiscard(allocator: std.mem.Allocator, atoms: []const model.Atom, bonds: []const model.Bond, graph: topology.Graph, rings: topology.RingMembership) !void {
    var result = try Fragmentation.init(allocator, atoms, bonds, graph, rings);
    result.deinit();
}

test "acyclic rotatable bonds split stable fragments and assign a parent tree" {
    const atoms = [_]model.Atom{ testAtom(0), testAtom(1), testAtom(2), testAtom(3), testAtom(4) };
    const bonds = [_]model.Bond{
        testBond(0, 0, 1),
        testBond(1, 1, 2),
        testBond(2, 2, 3),
        testBond(3, 3, 4),
    };
    var graph = try topology.Graph.init(std.testing.allocator, &atoms, &bonds);
    defer graph.deinit();
    var rings = try topology.RingMembership.init(std.testing.allocator, graph, &bonds);
    defer rings.deinit();
    var result = try Fragmentation.init(std.testing.allocator, &atoms, &bonds, graph, rings);
    defer result.deinit();

    try std.testing.expectEqual(@as(usize, 3), result.fragments.len);
    try std.testing.expectEqualSlices(core.ids.AtomId, &.{ core.ids.AtomId.fromIndex(0), core.ids.AtomId.fromIndex(1) }, result.members(core.ids.FragmentId.fromIndex(0)));
    try std.testing.expectEqualSlices(core.ids.AtomId, &.{core.ids.AtomId.fromIndex(2)}, result.members(core.ids.FragmentId.fromIndex(1)));
    try std.testing.expectEqualSlices(core.ids.AtomId, &.{ core.ids.AtomId.fromIndex(3), core.ids.AtomId.fromIndex(4) }, result.members(core.ids.FragmentId.fromIndex(2)));
    try std.testing.expectEqual(core.ids.FragmentId.fromIndex(0), result.main_fragments[0]);
    try std.testing.expectEqual(core.ids.FragmentId.fromIndex(0), result.fragments[1].parent);
    try std.testing.expectEqual(core.ids.FragmentId.fromIndex(1), result.fragments[2].parent);
    try std.testing.expect(result.fragments[0].flags.chain);
    try core.oom.checkAllocationFailures(std.testing.allocator, fragmentAndDiscard, .{ &atoms, &bonds, graph, rings });
}

test "ordinary rings remain one rigid fragment" {
    const atoms = [_]model.Atom{ testAtom(0), testAtom(1), testAtom(2), testAtom(3), testAtom(4), testAtom(5) };
    const bonds = [_]model.Bond{
        testBond(0, 0, 1),
        testBond(1, 1, 2),
        testBond(2, 2, 3),
        testBond(3, 3, 4),
        testBond(4, 4, 5),
        testBond(5, 5, 0),
    };
    var graph = try topology.Graph.init(std.testing.allocator, &atoms, &bonds);
    defer graph.deinit();
    var rings = try topology.RingMembership.init(std.testing.allocator, graph, &bonds);
    defer rings.deinit();
    var result = try Fragmentation.init(std.testing.allocator, &atoms, &bonds, graph, rings);
    defer result.deinit();
    try std.testing.expectEqual(@as(usize, 1), result.fragments.len);
    try std.testing.expectEqual(@as(u32, 1), result.fragments[0].ring_count);
    try std.testing.expect(!result.fragments[0].flags.chain);
}

test "prepared path fragments retain pinned canonical atom and creation order" {
    const input_atoms = [_]topology.prepare.TestAtom{ .{}, .{}, .{}, .{}, .{} };
    const input_bonds = [_]topology.prepare.TestBond{
        .{ .start = 0, .end = 1 },
        .{ .start = 1, .end = 2 },
        .{ .start = 2, .end = 3 },
        .{ .start = 3, .end = 4 },
    };
    var prepared = try topology.prepareInput(std.testing.allocator, topology.prepare.TestInput{ .atoms = &input_atoms, .bonds = &input_bonds });
    defer prepared.deinit();
    var result = try Fragmentation.init(
        std.testing.allocator,
        prepared.working.atoms,
        prepared.working.bonds,
        prepared.graph,
        prepared.rings,
    );
    defer result.deinit();
    try std.testing.expectEqualSlices(u32, &.{ 2, 1, 3, 0, 4 }, prepared.working.order.activeInternalToInput(5));
    try std.testing.expectEqual(@as(usize, 3), result.fragments.len);
    try std.testing.expectEqualSlices(core.ids.AtomId, &.{ core.ids.AtomId.fromIndex(1), core.ids.AtomId.fromIndex(3) }, result.members(core.ids.FragmentId.fromIndex(0)));
    try std.testing.expectEqualSlices(core.ids.AtomId, &.{core.ids.AtomId.fromIndex(0)}, result.members(core.ids.FragmentId.fromIndex(1)));
    try std.testing.expectEqualSlices(core.ids.AtomId, &.{ core.ids.AtomId.fromIndex(2), core.ids.AtomId.fromIndex(4) }, result.members(core.ids.FragmentId.fromIndex(2)));
    for (result.fragments) |fragment| try std.testing.expect(fragment.flags.chain);
    try std.testing.expectEqual(core.ids.FragmentId.fromIndex(0), result.fragments[1].parent);
    try std.testing.expectEqual(core.ids.FragmentId.fromIndex(1), result.fragments[2].parent);
}
