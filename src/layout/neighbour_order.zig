const std = @import("std");
const core = @import("core");
const model = @import("model");
const topology = @import("topology");

/// Upstream takes neighbours in input order except at degree four, where
/// `sketcherMinimizerAtom::orderAtomPriorities` ranks them. Coordinate layout
/// and DOF construction must share this order because upstream performs both
/// during the same traversal.
pub fn orderNeighbours(
    allocator: std.mem.Allocator,
    atoms: []const model.Atom,
    bonds: []const model.Bond,
    graph: topology.Graph,
    membership: topology.RingMembership,
    center: core.ids.AtomId,
    out: []core.ids.AtomId,
) core.errors.Error!void {
    const neighbours = graph.neighbors(center);
    if (out.len != neighbours.len) return error.InvalidMapping;
    @memcpy(out, neighbours);
    if (neighbours.len != 4) return;

    var weights: [4]f32 = undefined;
    for (neighbours, 0..) |neighbor, index| {
        weights[index] = try neighbourPriorityWeight(allocator, atoms, bonds, graph, membership, center, neighbor);
    }

    var rest: [4]core.ids.AtomId = undefined;
    var rest_weights: [4]f32 = undefined;
    var rest_count: usize = 0;
    for (neighbours, 0..) |neighbor, index| {
        rest[rest_count] = neighbor;
        rest_weights[rest_count] = weights[index];
        rest_count += 1;
    }
    const first = takeLowestWeight(rest[0..rest_count], rest_weights[0..rest_count]);
    rest_count -= 1;
    const second = takeLowestWeight(rest[0..rest_count], rest_weights[0..rest_count]);
    rest_count -= 1;

    const center_element = atoms[center.index()].atomic_number;
    if (center_element != .sulfur and center_element != .phosphorus) {
        out[0] = first;
        out[1] = rest[0];
        out[2] = rest[1];
        out[3] = second;
    } else {
        out[0] = first;
        out[1] = rest[0];
        out[2] = second;
        out[3] = rest[1];
    }
}

fn takeLowestWeight(atoms: []core.ids.AtomId, weights: []f32) core.ids.AtomId {
    var lowest: usize = 0;
    for (weights, 0..) |weight, index| {
        if (weight < weights[lowest]) lowest = index;
    }
    const taken = atoms[lowest];
    var index = lowest;
    while (index + 1 < atoms.len) : (index += 1) {
        atoms[index] = atoms[index + 1];
        weights[index] = weights[index + 1];
    }
    return taken;
}

/// One upstream term is deliberately absent: -2000 for a neighbour marked
/// `isSharedAndInner` while the centre is not. This port does not model that
/// pointer-era fused-ring flag, so assigning it here would invent behavior.
fn neighbourPriorityWeight(
    allocator: std.mem.Allocator,
    atoms: []const model.Atom,
    bonds: []const model.Bond,
    graph: topology.Graph,
    membership: topology.RingMembership,
    center: core.ids.AtomId,
    neighbor: core.ids.AtomId,
) core.errors.Error!f32 {
    const branch = try graph.reachableExcluding(allocator, neighbor, center);
    defer allocator.free(branch);
    var weight: f32 = @floatFromInt(branch.len);

    if (bondBetween(graph, center, neighbor)) |id| {
        const order = bonds[id.index()].effective_order;
        if (order == .double) weight -= 0.25;
        if (atoms[center.index()].atomic_number == .sulfur and order == .double) weight += 2000;
        if (sharesRing(membership, center, neighbor)) weight += 500;
    }
    if (atoms[neighbor.index()].atomic_number == .carbon) weight += 0.5;
    if (atoms[neighbor.index()].atomic_number == .hydrogen) weight -= 0.5;
    if (atoms[neighbor.index()].stereo != .unspecified) weight += 10000;
    if (atoms[center.index()].cross_layout and graph.degree(neighbor) > 1) weight += 200;
    for (graph.incidentBonds(neighbor)) |incident| if (bonds[incident.index()].effective_order == .double) {
        weight += 100;
        break;
    };
    return weight;
}

fn bondBetween(graph: topology.Graph, atom: core.ids.AtomId, other: core.ids.AtomId) ?core.ids.BondId {
    for (graph.neighbors(atom), graph.incidentBonds(atom)) |neighbor, bond| {
        if (neighbor == other) return bond;
    }
    return null;
}

fn sharesRing(membership: topology.RingMembership, first: core.ids.AtomId, second: core.ids.AtomId) bool {
    for (membership.atomRings(first)) |first_ring| {
        for (membership.atomRings(second)) |second_ring| if (first_ring == second_ring) return true;
    }
    return false;
}
