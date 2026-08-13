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
        try placeFragmentRings(allocator, atoms, bonds, membership, fragmentation, fragment, placed);

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
    try assembleFragments(atoms, bonds, graph, fragmentation);
    alignConstrainedMainFragments(atoms, fragmentation);
    fallbackOn3dCoordinates(atoms, fragmentation);
    restoreFixedCoordinates(atoms);
}

fn placeFragmentRings(
    allocator: std.mem.Allocator,
    atoms: []model.Atom,
    bonds: []const model.Bond,
    membership: topology.RingMembership,
    fragmentation: fragments.Fragmentation,
    fragment: fragments.Fragment,
    placed: []bool,
) core.errors.Error!void {
    var remaining: usize = fragment.ring_count;
    while (remaining != 0) {
        var selected: ?model.Ring = null;
        var selected_shared: usize = 0;
        var selected_score: usize = 0;
        for (membership.rings) |ring| {
            const ring_atoms = membership.atoms(ring.id);
            if (ring_atoms.len == 0 or ring_atoms.len >= 9 or
                fragmentation.atom_fragment[ring_atoms[0].index()] != fragment.id) continue;
            var already_complete = true;
            var shared: usize = 0;
            for (ring_atoms) |atom| {
                already_complete = already_complete and placed[atom.index()];
                shared += @intFromBool(placed[atom.index()]);
            }
            if (already_complete) continue;
            const score = shared * 10000 + @intFromBool(ring_atoms.len == 6) * 10 + ring_atoms.len;
            if (selected == null or score > selected_score) {
                selected = ring;
                selected_shared = shared;
                selected_score = score;
            }
        }
        const ring = selected orelse break;
        const ordered = try orderRing(allocator, ring, bonds, membership);
        defer allocator.free(ordered);
        const local = try regularRingCoordinates(allocator, ordered.len);
        defer allocator.free(local);
        if (selected_shared >= 2) {
            try alignFusedRing(atoms, ordered, local, placed);
        } else if (selected_shared == 1) {
            var pivot: usize = 0;
            for (ordered, 0..) |atom, index| if (placed[atom.index()]) {
                pivot = index;
                break;
            };
            const source_center = coordinateCenter(local);
            const existing_center = placedCenter(atoms, placed);
            const target = atoms[ordered[pivot].index()].coordinates;
            var outward = subtract(target, existing_center);
            if (length(outward) < 0.0001) outward = .{ .x = 1 };
            var source_direction = subtract(source_center, local[pivot]);
            if (length(source_direction) < 0.0001) source_direction = .{ .x = 1 };
            const rotation = std.math.atan2(outward.y, outward.x) - std.math.atan2(source_direction.y, source_direction.x);
            for (ordered, local) |atom, coordinate| if (!placed[atom.index()]) {
                atoms[atom.index()].coordinates = transformFromPivot(coordinate, local[pivot], target, rotation);
                placed[atom.index()] = true;
            };
        } else {
            for (ordered, local) |atom, coordinate| {
                atoms[atom.index()].coordinates = coordinate;
                placed[atom.index()] = true;
            }
        }
        remaining -= 1;
    }
}

fn regularRingCoordinates(allocator: std.mem.Allocator, count: usize) core.errors.Error![]core.math.Vec2 {
    if (count < 3) return error.InvalidMapping;
    const result = allocator.alloc(core.math.Vec2, count) catch return error.OutOfMemory;
    errdefer allocator.free(result);
    var coordinate: core.math.Vec2 = .{};
    const step = 2 * std.math.pi / @as(f32, @floatFromInt(count));
    for (result, 0..) |*destination, index| {
        destination.* = coordinate;
        const angle = step * @as(f32, @floatFromInt(index));
        coordinate.x += @cos(angle) * bond_length;
        coordinate.y -= @sin(angle) * bond_length;
    }
    return result;
}

fn alignFusedRing(atoms: []model.Atom, ordered: []const core.ids.AtomId, local: []const core.math.Vec2, placed: []bool) core.errors.Error!void {
    var first: ?usize = null;
    var last: usize = 0;
    for (ordered, 0..) |atom, index| if (placed[atom.index()]) {
        if (first == null) first = index;
        last = index;
    };
    const first_index = first orelse return error.InvalidMapping;
    if (first_index == last) return error.InvalidMapping;
    const target_first = atoms[ordered[first_index].index()].coordinates;
    const target_last = atoms[ordered[last].index()].coordinates;
    const source_first = local[first_index];
    const source_last = local[last];
    const target_angle = std.math.atan2(target_last.y - target_first.y, target_last.x - target_first.x);
    const source_angle = std.math.atan2(source_last.y - source_first.y, source_last.x - source_first.x);
    const rotation = target_angle - source_angle;
    var first_score: f32 = 0;
    var mirror_score: f32 = 0;
    const existing_center = placedCenter(atoms, placed);
    for (ordered, local) |atom, coordinate| {
        if (placed[atom.index()]) continue;
        const candidate = transformFromPivot(coordinate, source_first, target_first, rotation);
        const mirror = reflectAcrossLine(candidate, target_first, target_last);
        first_score += distance(candidate, existing_center);
        mirror_score += distance(mirror, existing_center);
    }
    const use_mirror = mirror_score > first_score;
    for (ordered, local) |atom, coordinate| if (!placed[atom.index()]) {
        const candidate = transformFromPivot(coordinate, source_first, target_first, rotation);
        atoms[atom.index()].coordinates = if (use_mirror) reflectAcrossLine(candidate, target_first, target_last) else candidate;
        placed[atom.index()] = true;
    };
}

fn transformFromPivot(point: core.math.Vec2, source: core.math.Vec2, target: core.math.Vec2, rotation: f32) core.math.Vec2 {
    const local = subtract(point, source);
    const cosine = @cos(rotation);
    const sine = @sin(rotation);
    return .{
        .x = target.x + local.x * cosine - local.y * sine,
        .y = target.y + local.x * sine + local.y * cosine,
    };
}

fn reflectAcrossLine(point: core.math.Vec2, line_start: core.math.Vec2, line_end: core.math.Vec2) core.math.Vec2 {
    const direction = subtract(line_end, line_start);
    const denominator = direction.x * direction.x + direction.y * direction.y;
    if (denominator < 0.000001) return point;
    const relative = subtract(point, line_start);
    const projection = (relative.x * direction.x + relative.y * direction.y) / denominator;
    const on_line = add(line_start, scale(direction, projection));
    return subtract(scale(on_line, 2), point);
}

fn placedCenter(atoms: []const model.Atom, placed: []const bool) core.math.Vec2 {
    var center: core.math.Vec2 = .{};
    var count: usize = 0;
    for (atoms, placed) |atom, is_placed| if (is_placed) {
        center = add(center, atom.coordinates);
        count += 1;
    };
    return if (count == 0) center else scale(center, 1 / @as(f32, @floatFromInt(count)));
}

fn coordinateCenter(coordinates: []const core.math.Vec2) core.math.Vec2 {
    var center: core.math.Vec2 = .{};
    for (coordinates) |coordinate| center = add(center, coordinate);
    return scale(center, 1 / @as(f32, @floatFromInt(coordinates.len)));
}

fn assembleFragments(atoms: []model.Atom, bonds: []const model.Bond, graph: topology.Graph, fragmentation: fragments.Fragmentation) core.errors.Error!void {
    _ = graph;
    for (1..fragmentation.fragments.len) |depth| {
        for (fragmentation.fragments) |fragment| {
            if (!fragment.parent.isValid() or fragmentDepth(fragmentation.fragments, fragment.id) != depth) continue;
            if (!fragment.bond_to_parent.isValid()) return error.InvalidMapping;
            const bond = bonds[fragment.bond_to_parent.index()];
            const child_endpoint = if (fragmentation.atom_fragment[bond.start.index()] == fragment.id) bond.start else bond.end;
            const parent_endpoint = if (child_endpoint == bond.start) bond.end else bond.start;
            const child_center = fragmentCenter(atoms, fragmentation.members(fragment.id));
            const parent_center = fragmentCenter(atoms, fragmentation.members(fragment.parent));
            var outward = subtract(atoms[parent_endpoint.index()].coordinates, parent_center);
            if (length(outward) < 0.0001) outward = .{ .x = 1 };
            outward = scale(outward, bond_length / length(outward));
            const target = add(atoms[parent_endpoint.index()].coordinates, outward);
            var child_direction = subtract(child_center, atoms[child_endpoint.index()].coordinates);
            if (length(child_direction) < 0.0001) child_direction = .{ .x = 1 };
            const angle = std.math.atan2(outward.y, outward.x) - std.math.atan2(child_direction.y, child_direction.x);
            const source = atoms[child_endpoint.index()].coordinates;
            for (fragmentation.members(fragment.id)) |atom| {
                atoms[atom.index()].coordinates = transformFromPivot(atoms[atom.index()].coordinates, source, target, angle);
            }
        }
    }
}

fn fragmentDepth(records: []const fragments.Fragment, fragment: core.ids.FragmentId) usize {
    var depth: usize = 0;
    var cursor = fragment;
    while (records[cursor.index()].parent.isValid()) {
        depth += 1;
        cursor = records[cursor.index()].parent;
        std.debug.assert(depth < records.len);
    }
    return depth;
}

fn alignConstrainedMainFragments(atoms: []model.Atom, fragmentation: fragments.Fragmentation) void {
    for (fragmentation.main_fragments) |main| {
        if (!main.isValid()) continue;
        const fragment = fragmentation.fragments[main.index()];
        if (!fragment.flags.constrained or fragment.flags.fixed) continue;
        const members = fragmentation.members(main);
        var source_center: core.math.Vec2 = .{};
        var target_center: core.math.Vec2 = .{};
        var constrained_count: usize = 0;
        for (members) |atom| if (atoms[atom.index()].constrained) {
            const target = atoms[atom.index()].template_coordinates orelse continue;
            source_center = add(source_center, atoms[atom.index()].coordinates);
            target_center = add(target_center, target);
            constrained_count += 1;
        };
        if (constrained_count == 0) continue;
        source_center = scale(source_center, 1 / @as(f32, @floatFromInt(constrained_count)));
        target_center = scale(target_center, 1 / @as(f32, @floatFromInt(constrained_count)));
        var dot: f32 = 0;
        var cross: f32 = 0;
        for (members) |atom| if (atoms[atom.index()].constrained) {
            const target = atoms[atom.index()].template_coordinates orelse continue;
            const source_delta = subtract(atoms[atom.index()].coordinates, source_center);
            const target_delta = subtract(target, target_center);
            dot += source_delta.x * target_delta.x + source_delta.y * target_delta.y;
            cross += source_delta.x * target_delta.y - source_delta.y * target_delta.x;
        };
        const angle = if (constrained_count > 1) std.math.atan2(cross, dot) else 0;
        for (members) |atom| {
            atoms[atom.index()].coordinates = transformFromPivot(atoms[atom.index()].coordinates, source_center, target_center, angle);
        }
    }
}

fn restoreFixedCoordinates(atoms: []model.Atom) void {
    for (atoms) |*atom| if (atom.fixed) {
        if (atom.template_coordinates) |coordinate| atom.coordinates = coordinate;
    };
}

fn fallbackOn3dCoordinates(atoms: []model.Atom, fragmentation: fragments.Fragmentation) void {
    for (fragmentation.fragments) |fragment| {
        const members = fragmentation.members(fragment.id);
        var generated_finite = true;
        var source_valid = true;
        for (members) |atom| {
            generated_finite = generated_finite and atoms[atom.index()].coordinates.isFinite();
            source_valid = source_valid and atoms[atom.index()].coordinates_3d != null and atoms[atom.index()].coordinates_3d.?.isFinite();
        }
        if (generated_finite or !source_valid) continue;
        for (members) |atom| {
            const source = atoms[atom.index()].coordinates_3d.?;
            atoms[atom.index()].coordinates = .{ .x = source.x, .y = source.y };
        }
    }
}

fn fragmentCenter(atoms: []const model.Atom, members: []const core.ids.AtomId) core.math.Vec2 {
    var center: core.math.Vec2 = .{};
    for (members) |atom| center = add(center, atoms[atom.index()].coordinates);
    return scale(center, 1 / @as(f32, @floatFromInt(members.len)));
}

fn add(left: core.math.Vec2, right: core.math.Vec2) core.math.Vec2 {
    return .{ .x = left.x + right.x, .y = left.y + right.y };
}

fn subtract(left: core.math.Vec2, right: core.math.Vec2) core.math.Vec2 {
    return .{ .x = left.x - right.x, .y = left.y - right.y };
}

fn scale(value: core.math.Vec2, factor: f32) core.math.Vec2 {
    return .{ .x = value.x * factor, .y = value.y * factor };
}

fn length(value: core.math.Vec2) f32 {
    return @sqrt(value.x * value.x + value.y * value.y);
}

fn distance(left: core.math.Vec2, right: core.math.Vec2) f32 {
    return length(subtract(left, right));
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

fn layoutFixture(atoms: []model.Atom, bonds: []const model.Bond) !void {
    var graph = try topology.Graph.init(std.testing.allocator, atoms, bonds);
    defer graph.deinit();
    var rings = try topology.RingMembership.init(std.testing.allocator, graph, bonds);
    defer rings.deinit();
    var split = try fragments.Fragmentation.init(std.testing.allocator, atoms, bonds, graph, rings);
    defer split.deinit();
    try initializeCoordinates(std.testing.allocator, atoms, bonds, graph, rings, split);
}

fn layoutAndDiscard(
    allocator: std.mem.Allocator,
    source_atoms: []const model.Atom,
    bonds: []const model.Bond,
    graph: topology.Graph,
    rings: topology.RingMembership,
    split: fragments.Fragmentation,
) !void {
    const atoms = try allocator.dupe(model.Atom, source_atoms);
    defer allocator.free(atoms);
    try initializeCoordinates(allocator, atoms, bonds, graph, rings, split);
}

fn expectBondLengths(atoms: []const model.Atom, bonds: []const model.Bond) !void {
    for (bonds) |bond| {
        try std.testing.expectApproxEqAbs(
            bond_length,
            distance(atoms[bond.start.index()].coordinates, atoms[bond.end.index()].coordinates),
            0.01,
        );
    }
}

test "fused rings align on their shared edge and extend outward" {
    var atoms: [6]model.Atom = undefined;
    for (&atoms, 0..) |*atom, index| atom.* = .{ .id = core.ids.AtomId.fromIndex(@intCast(index)), .input_index = @intCast(index), .atomic_number = .carbon };
    const pairs = [_][2]u32{ .{ 0, 1 }, .{ 1, 2 }, .{ 2, 3 }, .{ 3, 0 }, .{ 2, 5 }, .{ 5, 4 }, .{ 4, 3 } };
    var bonds: [pairs.len]model.Bond = undefined;
    for (&bonds, pairs, 0..) |*bond, pair, index| bond.* = .{
        .id = core.ids.BondId.fromIndex(@intCast(index)),
        .input_index = @intCast(index),
        .start = core.ids.AtomId.fromIndex(pair[0]),
        .end = core.ids.AtomId.fromIndex(pair[1]),
        .input_order = .single,
        .effective_order = .single,
    };
    try layoutFixture(&atoms, &bonds);
    try expectBondLengths(&atoms, &bonds);
    const first_center = scale(add(add(atoms[0].coordinates, atoms[1].coordinates), add(atoms[2].coordinates, atoms[3].coordinates)), 0.25);
    const second_center = scale(add(add(atoms[2].coordinates, atoms[3].coordinates), add(atoms[4].coordinates, atoms[5].coordinates)), 0.25);
    try std.testing.expect(distance(first_center, second_center) > bond_length * 0.5);

    var graph = try topology.Graph.init(std.testing.allocator, &atoms, &bonds);
    defer graph.deinit();
    var rings = try topology.RingMembership.init(std.testing.allocator, graph, &bonds);
    defer rings.deinit();
    var split = try fragments.Fragmentation.init(std.testing.allocator, &atoms, &bonds, graph, rings);
    defer split.deinit();
    try std.testing.checkAllAllocationFailures(
        std.testing.allocator,
        layoutAndDiscard,
        .{ &atoms, &bonds, graph, rings, split },
    );
}

test "spiro rings sharing one atom are placed on opposite sides" {
    var atoms: [5]model.Atom = undefined;
    for (&atoms, 0..) |*atom, index| atom.* = .{ .id = core.ids.AtomId.fromIndex(@intCast(index)), .input_index = @intCast(index), .atomic_number = .carbon };
    const pairs = [_][2]u32{ .{ 0, 1 }, .{ 1, 2 }, .{ 2, 0 }, .{ 0, 3 }, .{ 3, 4 }, .{ 4, 0 } };
    var bonds: [pairs.len]model.Bond = undefined;
    for (&bonds, pairs, 0..) |*bond, pair, index| bond.* = .{
        .id = core.ids.BondId.fromIndex(@intCast(index)),
        .input_index = @intCast(index),
        .start = core.ids.AtomId.fromIndex(pair[0]),
        .end = core.ids.AtomId.fromIndex(pair[1]),
        .input_order = .single,
        .effective_order = .single,
    };
    try layoutFixture(&atoms, &bonds);
    try expectBondLengths(&atoms, &bonds);
    const first_center = scale(add(add(atoms[0].coordinates, atoms[1].coordinates), atoms[2].coordinates), 1.0 / 3.0);
    const second_center = scale(add(add(atoms[0].coordinates, atoms[3].coordinates), atoms[4].coordinates), 1.0 / 3.0);
    const first_direction = subtract(first_center, atoms[0].coordinates);
    const second_direction = subtract(second_center, atoms[0].coordinates);
    try std.testing.expect(first_direction.x * second_direction.x + first_direction.y * second_direction.y < 0);
}

test "fragment assembly preserves every acyclic parent bond length" {
    var atoms: [5]model.Atom = undefined;
    for (&atoms, 0..) |*atom, index| atom.* = .{ .id = core.ids.AtomId.fromIndex(@intCast(index)), .input_index = @intCast(index), .atomic_number = .carbon };
    var bonds: [4]model.Bond = undefined;
    for (&bonds, 0..) |*bond, index| bond.* = .{
        .id = core.ids.BondId.fromIndex(@intCast(index)),
        .input_index = @intCast(index),
        .start = core.ids.AtomId.fromIndex(@intCast(index)),
        .end = core.ids.AtomId.fromIndex(@intCast(index + 1)),
        .input_order = .single,
        .effective_order = .single,
    };
    try layoutFixture(&atoms, &bonds);
    try expectBondLengths(&atoms, &bonds);
    for (atoms) |atom| try std.testing.expect(atom.coordinates.isFinite());
    const first = atoms;
    try layoutFixture(&atoms, &bonds);
    for (first, atoms) |left, right| try std.testing.expectEqual(left.coordinates, right.coordinates);
}

test "constrained alignment, fixed reset, and valid 3D fallback are deterministic" {
    var atoms: [3]model.Atom = undefined;
    for (&atoms, 0..) |*atom, index| atom.* = .{
        .id = core.ids.AtomId.fromIndex(@intCast(index)),
        .input_index = @intCast(index),
        .atomic_number = .carbon,
    };
    var bonds = [_]model.Bond{
        .{ .id = core.ids.BondId.fromIndex(0), .input_index = 0, .start = core.ids.AtomId.fromIndex(0), .end = core.ids.AtomId.fromIndex(1), .input_order = .single, .effective_order = .single },
        .{ .id = core.ids.BondId.fromIndex(1), .input_index = 1, .start = core.ids.AtomId.fromIndex(1), .end = core.ids.AtomId.fromIndex(2), .input_order = .single, .effective_order = .single },
    };
    atoms[0].constrained = true;
    atoms[0].template_coordinates = .{ .x = 10, .y = 20 };
    atoms[1].constrained = true;
    atoms[1].template_coordinates = .{ .x = 10, .y = 70 };
    try layoutFixture(&atoms, &bonds);
    try std.testing.expectApproxEqAbs(@as(f32, 10), atoms[0].coordinates.x, 0.001);
    try std.testing.expectApproxEqAbs(@as(f32, 20), atoms[0].coordinates.y, 0.001);
    try std.testing.expectApproxEqAbs(@as(f32, 10), atoms[1].coordinates.x, 0.001);
    try std.testing.expectApproxEqAbs(@as(f32, 70), atoms[1].coordinates.y, 0.001);

    atoms[2].fixed = true;
    atoms[2].template_coordinates = .{ .x = -12, .y = 34 };
    try layoutFixture(&atoms, &bonds);
    try std.testing.expectEqual(core.math.Vec2{ .x = -12, .y = 34 }, atoms[2].coordinates);

    for (&atoms, 0..) |*atom, index| {
        atom.fixed = false;
        atom.constrained = true;
        atom.template_coordinates = .{ .x = std.math.nan(f32), .y = 0 };
        atom.coordinates_3d = .{ .x = @floatFromInt(index * 3), .y = @floatFromInt(index * 5), .z = 7 };
    }
    try layoutFixture(&atoms, &bonds);
    for (atoms, 0..) |atom, index| {
        try std.testing.expectEqual(@as(f32, @floatFromInt(index * 3)), atom.coordinates.x);
        try std.testing.expectEqual(@as(f32, @floatFromInt(index * 5)), atom.coordinates.y);
    }
}
