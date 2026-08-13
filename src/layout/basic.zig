const std = @import("std");
const core = @import("core");
const model = @import("model");
const topology = @import("topology");
const fragments = @import("fragments.zig");

pub const bond_length: f32 = 50;

/// Generate deterministic local coordinates for ordinary rigid fragments.
/// Template and macrocycle dispatch are deliberately owned by their later
/// layers; this seam covers regular rings and the acyclic fallback.
pub fn initializeCoordinates(
    allocator: std.mem.Allocator,
    atoms: []model.Atom,
    bonds: []const model.Bond,
    graph: topology.Graph,
    membership: topology.RingMembership,
    fragmentation: fragments.Fragmentation,
) core.errors.Error!void {
    const placed = allocator.alloc(bool, atoms.len) catch return error.OutOfMemory;
    defer allocator.free(placed);
    @memset(placed, false);
    const queue = allocator.alloc(core.ids.AtomId, atoms.len) catch return error.OutOfMemory;
    defer allocator.free(queue);

    for (fragmentation.fragments) |fragment| {
        for (membership.rings) |ring| {
            const ring_atoms = membership.atoms(ring.id);
            if (ring_atoms.len == 0 or fragmentation.atom_fragment[ring_atoms[0].index()] != fragment.id) continue;
            if (ring_atoms.len >= 9) continue;
            const ordered = try orderRing(allocator, ring, bonds, membership);
            defer allocator.free(ordered);
            var coordinate: core.math.Vec2 = .{};
            const step = 2 * std.math.pi / @as(f32, @floatFromInt(ordered.len));
            for (ordered, 0..) |atom, index| {
                if (!placed[atom.index()]) atoms[atom.index()].coordinates = coordinate;
                placed[atom.index()] = true;
                const angle = step * @as(f32, @floatFromInt(index));
                coordinate.x += @cos(angle) * bond_length;
                coordinate.y -= @sin(angle) * bond_length;
            }
        }

        const members = fragmentation.members(fragment.id);
        var head: usize = 0;
        var tail: usize = 0;
        for (members) |atom| if (placed[atom.index()]) {
            queue[tail] = atom;
            tail += 1;
        };
        if (tail == 0) {
            const start = members[0];
            atoms[start.index()].coordinates = .{};
            placed[start.index()] = true;
            queue[0] = start;
            tail = 1;
        }
        while (head < tail) : (head += 1) {
            const center = queue[head];
            var unplaced_count: usize = 0;
            for (graph.neighbors(center)) |neighbor| {
                if (fragmentation.atom_fragment[neighbor.index()] == fragment.id and !placed[neighbor.index()]) unplaced_count += 1;
            }
            if (unplaced_count == 0) continue;
            const base_angle = parentAngle(atoms, graph, fragmentation, fragment.id, center);
            var generated: usize = 0;
            for (graph.neighbors(center)) |neighbor| {
                if (fragmentation.atom_fragment[neighbor.index()] != fragment.id or placed[neighbor.index()]) continue;
                const spread = if (unplaced_count == 1)
                    @as(f32, 0)
                else
                    (@as(f32, @floatFromInt(generated)) - @as(f32, @floatFromInt(unplaced_count - 1)) * 0.5) * (2 * std.math.pi / 3);
                const angle = base_angle + spread;
                atoms[neighbor.index()].coordinates = .{
                    .x = atoms[center.index()].coordinates.x + @cos(angle) * bond_length,
                    .y = atoms[center.index()].coordinates.y + @sin(angle) * bond_length,
                };
                placed[neighbor.index()] = true;
                queue[tail] = neighbor;
                tail += 1;
                generated += 1;
            }
        }
    }
}

fn orderRing(allocator: std.mem.Allocator, ring: model.Ring, bonds: []const model.Bond, membership: topology.RingMembership) core.errors.Error![]core.ids.AtomId {
    const members = membership.atoms(ring.id);
    const result = allocator.alloc(core.ids.AtomId, members.len) catch return error.OutOfMemory;
    errdefer allocator.free(result);
    result[0] = members[0];
    var previous = core.ids.AtomId.invalid;
    for (1..members.len) |index| {
        var next = core.ids.AtomId.invalid;
        for (membership.ringBonds(ring.id)) |bond_id| {
            const bond = bonds[bond_id.index()];
            const candidate = if (bond.start == result[index - 1]) bond.end else if (bond.end == result[index - 1]) bond.start else continue;
            if (candidate == previous) continue;
            if (index + 1 < members.len and std.mem.indexOfScalar(core.ids.AtomId, result[0..index], candidate) != null) continue;
            if (!next.isValid() or candidate.index() < next.index()) next = candidate;
        }
        if (!next.isValid()) return error.InvalidMapping;
        previous = result[index - 1];
        result[index] = next;
    }
    return result;
}

fn parentAngle(atoms: []const model.Atom, graph: topology.Graph, fragmentation: fragments.Fragmentation, fragment: core.ids.FragmentId, center: core.ids.AtomId) f32 {
    for (graph.neighbors(center)) |neighbor| {
        if (fragmentation.atom_fragment[neighbor.index()] != fragment or
            (atoms[neighbor.index()].coordinates.x == atoms[center.index()].coordinates.x and atoms[neighbor.index()].coordinates.y == atoms[center.index()].coordinates.y)) continue;
        const direction = atoms[center.index()].coordinates;
        const parent = atoms[neighbor.index()].coordinates;
        return std.math.atan2(direction.y - parent.y, direction.x - parent.x) + std.math.pi / 3;
    }
    return 0;
}

test "regular ring coordinate walk preserves the upstream bond length" {
    var atoms: [6]model.Atom = undefined;
    var bonds: [6]model.Bond = undefined;
    for (&atoms, 0..) |*atom, index| atom.* = .{ .id = core.ids.AtomId.fromIndex(@intCast(index)), .input_index = @intCast(index), .atomic_number = .carbon };
    for (&bonds, 0..) |*bond, index| bond.* = .{
        .id = core.ids.BondId.fromIndex(@intCast(index)),
        .input_index = @intCast(index),
        .start = core.ids.AtomId.fromIndex(@intCast(index)),
        .end = core.ids.AtomId.fromIndex(@intCast((index + 1) % atoms.len)),
        .input_order = .single,
        .effective_order = .single,
    };
    var graph = try topology.Graph.init(std.testing.allocator, &atoms, &bonds);
    defer graph.deinit();
    var rings = try topology.RingMembership.init(std.testing.allocator, graph, &bonds);
    defer rings.deinit();
    var split = try fragments.Fragmentation.init(std.testing.allocator, &atoms, &bonds, graph, rings);
    defer split.deinit();
    try initializeCoordinates(std.testing.allocator, &atoms, &bonds, graph, rings, split);
    for (bonds) |bond| {
        const a = atoms[bond.start.index()].coordinates;
        const b = atoms[bond.end.index()].coordinates;
        try std.testing.expectApproxEqAbs(bond_length, @sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y)), 0.001);
    }
}
