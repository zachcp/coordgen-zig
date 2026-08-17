//! Bend-interaction group construction (cgz-r33).
//!
//! `optimize.buildBendInteractions` consumes `BendGroup` values and has never
//! had a caller, because nothing built them. It lives here rather than in the
//! optimize layer because it needs the graph and ring structure, and the
//! approved module edges give `optimize` neither: `optimize` may import only
//! core, model and geometry. The generator layer is the first that can see
//! topology and optimize at once.
//!
//! Mirrors `CoordgenMinimizer::addBendInteractionsOfMolecule`. The angle
//! redistribution, the rigid-atom override and the four-neighbour cases are
//! already implemented by the consumer, so what is built here is exactly its
//! input: per centre, the clockwise-adjacent neighbour pairs with their
//! starting rest angle and ring context.

const std = @import("std");
const core = @import("core");
const model = @import("model");
const geometry = @import("geometry");
const topology = @import("topology");
const optimize = @import("optimize");

pub const Groups = struct {
    allocator: std.mem.Allocator,
    groups: []optimize.BendGroup,
    candidates: []optimize.BendCandidate,

    pub fn deinit(self: *Groups) void {
        self.allocator.free(self.candidates);
        self.allocator.free(self.groups);
        self.* = undefined;
    }
};

/// Build one group per atom with more than one neighbour, in atom order.
///
/// Candidate storage is one flat allocation; each group's `candidates` slice
/// borrows from it, so `Groups.deinit` frees both and the groups must not
/// outlive it.
pub fn build(
    allocator: std.mem.Allocator,
    atoms: []const model.Atom,
    bonds: []const model.Bond,
    graph: topology.Graph,
    membership: topology.RingMembership,
    analysis: topology.rings.Analysis,
) core.errors.Error!Groups {
    var groups: std.ArrayList(optimize.BendGroup) = .empty;
    defer groups.deinit(allocator);
    var candidates: std.ArrayList(optimize.BendCandidate) = .empty;
    defer candidates.deinit(allocator);
    // Ranges are recorded first and resolved to slices after the candidate
    // array has stopped growing, because growth invalidates earlier pointers.
    var ranges: std.ArrayList(struct { start: usize, len: usize }) = .empty;
    defer ranges.deinit(allocator);

    var ordered_buffer: [16]core.ids.AtomId = undefined;

    for (atoms) |atom| {
        const neighbors = graph.neighbors(atom.id);
        if (neighbors.len < 2) continue;
        const ordered = try topology.clockwiseOrderedNeighbors(
            graph,
            atoms,
            atom.id,
            &ordered_buffer,
        );

        const base_angle = try baseRestAngle(atom.id, ordered, bonds, graph);
        const start = candidates.items.len;
        var inverted_macrocycle_bond = false;

        for (ordered, 0..) |first, index| {
            // Two neighbours produce one interaction, not two: the i = 0 and
            // i = 1 pairs are the same pair.
            if (ordered.len == 2 and index == 1) continue;
            const second = ordered[(index + ordered.len - 1) % ordered.len];
            const ring = membership.sameRing(atom.id, first, second);
            var context: optimize.BendRingContext = .non_ring;
            if (ring) |ring_id| {
                if (membership.atoms(ring_id).len >= topology.rings.macrocycle_size) {
                    context = .macrocycle_ring;
                    if (ordered.len == 3 and try invertsMacrocycleBond(
                        atoms,
                        atom.id,
                        first,
                        second,
                        ordered,
                    )) {
                        inverted_macrocycle_bond = true;
                    }
                } else {
                    context = .{ .small_ring = effectiveRingSize(membership, analysis, ring_id) };
                }
            }
            candidates.append(allocator, .{
                .atom_a = first,
                .atom_b = second,
                .initial_rest_degrees = base_angle,
                .ring = context,
            }) catch return error.OutOfMemory;
        }

        groups.append(allocator, .{
            .center = atom.id,
            .candidates = &.{},
            .cross_layout = atom.cross_layout,
            .inverted_macrocycle_bond = inverted_macrocycle_bond,
        }) catch return error.OutOfMemory;
        ranges.append(allocator, .{
            .start = start,
            .len = candidates.items.len - start,
        }) catch return error.OutOfMemory;
    }

    const owned_candidates = candidates.toOwnedSlice(allocator) catch return error.OutOfMemory;
    errdefer allocator.free(owned_candidates);
    const owned_groups = groups.toOwnedSlice(allocator) catch return error.OutOfMemory;
    for (owned_groups, ranges.items) |*group, range| {
        group.candidates = owned_candidates[range.start .. range.start + range.len];
    }
    return .{
        .allocator = allocator,
        .groups = owned_groups,
        .candidates = owned_candidates,
    };
}

/// 120 degrees for two neighbours, 180 when their two bond orders sum above
/// three, and an even split otherwise.
fn baseRestAngle(
    center: core.ids.AtomId,
    ordered: []const core.ids.AtomId,
    bonds: []const model.Bond,
    graph: topology.Graph,
) core.errors.Error!f32 {
    if (ordered.len > 2) return 360 / @as(f32, @floatFromInt(ordered.len));
    var order_sum: u32 = 0;
    for (graph.incidentBonds(center)) |bond_id| {
        if (bond_id.index() >= bonds.len) return error.InvalidMapping;
        order_sum += @backingInt(bonds[bond_id.index()].effective_order);
    }
    return if (order_sum > 3) 180 else 120;
}

/// Upstream's three-neighbour macrocycle test: the bend inverts when the far
/// neighbour lies on the same side of the centre-to-first line as the other
/// neighbour of the pair.
fn invertsMacrocycleBond(
    atoms: []const model.Atom,
    center: core.ids.AtomId,
    first: core.ids.AtomId,
    second: core.ids.AtomId,
    ordered: []const core.ids.AtomId,
) core.errors.Error!bool {
    var other: ?core.ids.AtomId = null;
    for (ordered) |candidate| {
        if (candidate == first or candidate == second) continue;
        other = candidate;
        break;
    }
    const third = other orelse return false;
    if (center.index() >= atoms.len or first.index() >= atoms.len or
        second.index() >= atoms.len or third.index() >= atoms.len)
    {
        return error.InvalidMapping;
    }
    return geometry.sameSide(
        atoms[second.index()].coordinates,
        atoms[third.index()].coordinates,
        atoms[first.index()].coordinates,
        atoms[center.index()].coordinates,
    );
}

/// Rings drawn as fused occupy a bigger effective ring, so the rest angle is
/// computed against the inflated size: for every non-macrocyclic ring fused
/// to this one across more than two atoms, add the atoms that fusion brings.
fn effectiveRingSize(
    membership: topology.RingMembership,
    analysis: topology.rings.Analysis,
    ring: core.ids.RingId,
) u32 {
    var size: u32 = @intCast(membership.atoms(ring).len);
    for (analysis.fusedWith(ring)) |fusion| {
        const other_atoms = membership.atoms(fusion.other);
        if (other_atoms.len >= topology.rings.macrocycle_size) continue;
        const shared = analysis.fusionAtoms(fusion).len;
        if (shared <= 2) continue;
        size += @intCast(other_atoms.len - shared);
    }
    return size;
}

const testing = std.testing;

fn testBond(index: u32, start: u32, end: u32, order: core.chemistry.BondOrder) model.Bond {
    return .{
        .id = core.ids.BondId.fromIndex(index),
        .input_index = index,
        .start = core.ids.AtomId.fromIndex(start),
        .end = core.ids.AtomId.fromIndex(end),
        .input_order = order,
        .effective_order = order,
    };
}

fn testAtom(index: u32, x: f32, y: f32) model.Atom {
    return .{
        .id = core.ids.AtomId.fromIndex(index),
        .input_index = index,
        .atomic_number = .carbon,
        .coordinates = .{ .x = x, .y = y },
    };
}

test "a two-neighbour centre emits one candidate, and bond order decides 120 against 180" {
    const atoms = [_]model.Atom{
        testAtom(0, 0, 0),
        testAtom(1, -50, 0),
        testAtom(2, 50, 0),
    };
    const single = [_]model.Bond{
        testBond(0, 0, 1, .single),
        testBond(1, 0, 2, .single),
    };
    var graph = try topology.Graph.init(testing.allocator, &atoms, &single);
    defer graph.deinit();
    var membership = try topology.RingMembership.init(testing.allocator, graph, &single);
    defer membership.deinit();
    var analysis = try topology.rings.Analysis.init(testing.allocator, membership, &atoms, &single);
    defer analysis.deinit();

    var groups = try build(testing.allocator, &atoms, &single, graph, membership, analysis);
    defer groups.deinit();
    // Only the centre has two neighbours.
    try testing.expectEqual(@as(usize, 1), groups.groups.len);
    try testing.expectEqual(core.ids.AtomId.fromIndex(0), groups.groups[0].center);
    // One interaction, not two: the first and last pair coincide.
    try testing.expectEqual(@as(usize, 1), groups.groups[0].candidates.len);
    try testing.expectEqual(@as(f32, 120), groups.groups[0].candidates[0].initial_rest_degrees);
    try testing.expect(groups.groups[0].candidates[0].ring == .non_ring);

    // A double plus a triple sums above three and straightens the centre.
    const unsaturated = [_]model.Bond{
        testBond(0, 0, 1, .double),
        testBond(1, 0, 2, .triple),
    };
    var graph2 = try topology.Graph.init(testing.allocator, &atoms, &unsaturated);
    defer graph2.deinit();
    var membership2 = try topology.RingMembership.init(testing.allocator, graph2, &unsaturated);
    defer membership2.deinit();
    var analysis2 = try topology.rings.Analysis.init(testing.allocator, membership2, &atoms, &unsaturated);
    defer analysis2.deinit();
    var straight = try build(testing.allocator, &atoms, &unsaturated, graph2, membership2, analysis2);
    defer straight.deinit();
    try testing.expectEqual(@as(f32, 180), straight.groups[0].candidates[0].initial_rest_degrees);
}

test "a ring centre reports its ring context and an even split for three neighbours" {
    // Cyclopropane with one exocyclic substituent on atom 0.
    const atoms = [_]model.Atom{
        testAtom(0, 0, 0),
        testAtom(1, 50, 0),
        testAtom(2, 25, 43),
        testAtom(3, -50, 0),
    };
    const bonds = [_]model.Bond{
        testBond(0, 0, 1, .single),
        testBond(1, 1, 2, .single),
        testBond(2, 2, 0, .single),
        testBond(3, 0, 3, .single),
    };
    var graph = try topology.Graph.init(testing.allocator, &atoms, &bonds);
    defer graph.deinit();
    var membership = try topology.RingMembership.init(testing.allocator, graph, &bonds);
    defer membership.deinit();
    var analysis = try topology.rings.Analysis.init(testing.allocator, membership, &atoms, &bonds);
    defer analysis.deinit();

    var groups = try build(testing.allocator, &atoms, &bonds, graph, membership, analysis);
    defer groups.deinit();

    var center_group: ?optimize.BendGroup = null;
    for (groups.groups) |group| {
        if (group.center == core.ids.AtomId.fromIndex(0)) center_group = group;
    }
    const group = center_group orelse return error.TestUnexpectedResult;
    // Three neighbours: three candidates at 360/3.
    try testing.expectEqual(@as(usize, 3), group.candidates.len);
    for (group.candidates) |candidate| {
        try testing.expectEqual(@as(f32, 120), candidate.initial_rest_degrees);
    }
    // Exactly one of the three pairs is the in-ring pair 1-0-2.
    var ring_pairs: usize = 0;
    for (group.candidates) |candidate| {
        switch (candidate.ring) {
            .small_ring => |size| {
                ring_pairs += 1;
                // Unfused, so the effective size is the ring's own size.
                try testing.expectEqual(@as(u32, 3), size);
            },
            else => {},
        }
    }
    try testing.expectEqual(@as(usize, 1), ring_pairs);
}
