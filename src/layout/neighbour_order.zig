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
    analysis: topology.rings.Analysis,
    center: core.ids.AtomId,
    out: []core.ids.AtomId,
) core.errors.Error!void {
    const neighbours = graph.neighbors(center);
    if (out.len != neighbours.len) return error.InvalidMapping;
    @memcpy(out, neighbours);
    if (neighbours.len != 4) return;

    var weights: [4]f32 = undefined;
    for (neighbours, 0..) |neighbor, index| {
        weights[index] = try topology.stereo.atomPriorityWeight(
            allocator,
            atoms,
            bonds,
            graph,
            membership,
            analysis.shared_and_inner,
            center,
            neighbor,
        );
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
