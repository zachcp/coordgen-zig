const std = @import("std");
const core = @import("core");
const geometry = @import("geometry");

pub const PenaltyContext = struct {
    constrained_flip: bool = false,
    chain_with_chain_parent: bool = false,
};

pub fn affectedAtoms(collection: core.dof.Collection, dof: core.dof.Dof) core.errors.Error![]const core.ids.AtomId {
    const start = dof.affected_atoms.start;
    const end = std.math.add(u32, start, dof.affected_atoms.len) catch return error.InvalidMapping;
    if (end > collection.affected_atoms.len) return error.InvalidMapping;
    return collection.affected_atoms[start..end];
}

/// Apply one DOF to fragment-local coordinates. Callers rebuild global
/// fragment placement after applying all current states, as pinned upstream
/// does in CoordgenDOFSolutions::scoreCurrentSolution.
pub fn apply(
    dof: core.dof.Dof,
    affected: []const core.ids.AtomId,
    coordinates: []core.math.Vec2,
) !void {
    try dof.state.validate();
    for (affected) |atom| if (atom.index() >= coordinates.len) return error.InvalidAtomIndex;
    if (dof.state.current == 0) return;
    switch (dof.payload) {
        .flip_fragment => {
            for (affected) |atom| coordinates[atom.index()].y = -coordinates[atom.index()].y;
        },
        .change_parent_bond_length => {
            const factor = stateFactor(1.6, dof.state.current);
            const move_by = core.math.bond_length * (factor - 1);
            for (affected) |atom| coordinates[atom.index()].x += move_by;
        },
        .rotate_fragment => {
            var angle = @as(f32, std.math.pi) / 180 * 15 * @as(f32, @floatFromInt((dof.state.current + 1) / 2));
            if (dof.state.current % 2 == 0) angle = -angle;
            const origin = core.math.Vec2{ .x = -core.math.bond_length };
            for (affected) |atom| {
                const offset = subtract(coordinates[atom.index()], origin);
                coordinates[atom.index()] = add(geometry.rotate(offset, @sin(angle), @cos(angle)), origin);
            }
        },
        .scale_atoms => |payload| {
            if (payload.pivot.index() >= coordinates.len) return error.InvalidAtomIndex;
            const pivot = coordinates[payload.pivot.index()];
            for (affected) |atom| coordinates[atom.index()] = add(pivot, scale(subtract(coordinates[atom.index()], pivot), 0.4));
        },
        .scale_fragment => {
            const factor = stateFactor(1.4, dof.state.current);
            for (affected) |atom| coordinates[atom.index()] = scale(coordinates[atom.index()], factor);
        },
        .invert_bond => |payload| {
            if (payload.pivot.index() >= coordinates.len or payload.bound.index() >= coordinates.len) return error.InvalidAtomIndex;
            const pivot = coordinates[payload.pivot.index()];
            const bond = subtract(coordinates[payload.bound.index()], pivot);
            const line_a = add(pivot, .{ .x = bond.y, .y = -bond.x });
            const line_b = subtract(pivot, .{ .x = bond.y, .y = -bond.x });
            for (affected) |atom| coordinates[atom.index()] = reflect(coordinates[atom.index()], line_a, line_b);
        },
        .flip_ring => |payload| {
            if (payload.pivot_a.index() >= coordinates.len or payload.pivot_b.index() >= coordinates.len) return error.InvalidAtomIndex;
            const line_a = coordinates[payload.pivot_a.index()];
            const line_b = coordinates[payload.pivot_b.index()];
            for (affected) |atom| coordinates[atom.index()] = reflect(coordinates[atom.index()], line_a, line_b);
        },
    }
}

pub fn penalty(dof: core.dof.Dof, affected_count: usize, context: PenaltyContext) !f32 {
    try dof.state.validate();
    const state_level = @as(f32, @floatFromInt((dof.state.current + 1) / 2));
    return switch (dof.payload) {
        .flip_fragment => @as(f32, @floatFromInt(@as(u32, @intFromBool(context.chain_with_chain_parent)) * 10 +
            @as(u32, @intFromBool(dof.state.current != 0 and context.constrained_flip)) * 1000)),
        .change_parent_bond_length => if (dof.state.current == 0) 0 else 200 * state_level,
        .rotate_fragment => if (dof.state.current == 0) 0 else 400 * state_level,
        .scale_atoms => if (dof.state.current == 0) 0 else 50 * @as(f32, @floatFromInt(affected_count)),
        .scale_fragment => if (dof.state.current == 0) 0 else 500 * state_level,
        .invert_bond => if (dof.state.current == 0) 0 else 100,
        .flip_ring => |payload| if (dof.state.current == 0) 0 else 200 * @as(f32, @floatFromInt(payload.penalty_multiplier)),
    };
}

fn stateFactor(base: f32, state: u32) f32 {
    var factor = std.math.pow(f32, base, @floatFromInt((state + 1) / 2));
    if (state % 2 == 0) factor = 1 / factor;
    return factor;
}

fn reflect(point: core.math.Vec2, line_a: core.math.Vec2, line_b: core.math.Vec2) core.math.Vec2 {
    const direction = subtract(line_b, line_a);
    const denominator = direction.x * direction.x + direction.y * direction.y;
    if (denominator == 0) return point;
    const relative = subtract(point, line_a);
    const projection = (relative.x * direction.x + relative.y * direction.y) / denominator;
    const on_line = add(line_a, scale(direction, projection));
    return subtract(scale(on_line, 2), point);
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

fn testDof(payload: core.dof.Payload, state: u32, count: u32) core.dof.Dof {
    return .{
        .id = core.ids.DofId.fromIndex(0),
        .fragment = core.ids.FragmentId.fromIndex(0),
        .state = .{ .current = state, .optimal = 0, .count = count, .tier = 0 },
        .payload = payload,
    };
}

test "all seven DOF applications preserve pinned local transforms" {
    const affected = [_]core.ids.AtomId{core.ids.AtomId.fromIndex(2)};

    var coordinates = [_]core.math.Vec2{ .{}, .{}, .{ .x = 2, .y = 3 } };
    try apply(testDof(.{ .flip_fragment = .{} }, 1, 2), &affected, &coordinates);
    try std.testing.expectEqual(core.math.Vec2{ .x = 2, .y = -3 }, coordinates[2]);

    coordinates[2] = .{ .x = 2, .y = 3 };
    try apply(testDof(.{ .scale_fragment = .{} }, 1, 5), &affected, &coordinates);
    try std.testing.expectApproxEqAbs(@as(f32, 2.8), coordinates[2].x, 0.0001);

    coordinates[1] = .{ .x = 1, .y = 1 };
    coordinates[2] = .{ .x = 6, .y = 1 };
    try apply(testDof(.{ .scale_atoms = .{ .pivot = core.ids.AtomId.fromIndex(1) } }, 1, 2), &affected, &coordinates);
    try std.testing.expectEqual(core.math.Vec2{ .x = 3, .y = 1 }, coordinates[2]);

    coordinates[2] = .{ .x = 0, .y = 0 };
    try apply(testDof(.{ .change_parent_bond_length = .{} }, 1, 7), &affected, &coordinates);
    try std.testing.expectApproxEqAbs(@as(f32, 30), coordinates[2].x, 0.0001);

    coordinates[2] = .{ .x = 0, .y = 0 };
    try apply(testDof(.{ .rotate_fragment = .{} }, 1, 5), &affected, &coordinates);
    try std.testing.expectApproxEqAbs(@as(f32, -1.7037087), coordinates[2].x, 0.0001);
    try std.testing.expectApproxEqAbs(@as(f32, -12.940952), coordinates[2].y, 0.0001);

    coordinates[0] = .{};
    coordinates[1] = .{ .x = 1 };
    coordinates[2] = .{ .x = 2, .y = 3 };
    try apply(testDof(.{ .invert_bond = .{ .pivot = core.ids.AtomId.fromIndex(0), .bound = core.ids.AtomId.fromIndex(1) } }, 1, 2), &affected, &coordinates);
    try std.testing.expectEqual(core.math.Vec2{ .x = -2, .y = 3 }, coordinates[2]);

    coordinates[0] = .{};
    coordinates[1] = .{ .x = 1 };
    coordinates[2] = .{ .x = 2, .y = 3 };
    try apply(testDof(.{ .flip_ring = .{
        .ring = core.ids.RingId.fromIndex(0),
        .pivot_a = core.ids.AtomId.fromIndex(0),
        .pivot_b = core.ids.AtomId.fromIndex(1),
        .penalty_multiplier = 3,
    } }, 1, 2), &affected, &coordinates);
    try std.testing.expectEqual(core.math.Vec2{ .x = 2, .y = -3 }, coordinates[2]);
}

test "DOF penalties preserve state levels and flip context" {
    try std.testing.expectEqual(@as(f32, 1010), try penalty(testDof(.{ .flip_fragment = .{} }, 1, 2), 0, .{
        .constrained_flip = true,
        .chain_with_chain_parent = true,
    }));
    try std.testing.expectEqual(@as(f32, 400), try penalty(testDof(.{ .change_parent_bond_length = .{} }, 3, 7), 0, .{}));
    try std.testing.expectEqual(@as(f32, 800), try penalty(testDof(.{ .rotate_fragment = .{} }, 3, 5), 0, .{}));
    try std.testing.expectEqual(@as(f32, 150), try penalty(testDof(.{ .scale_atoms = .{ .pivot = core.ids.AtomId.fromIndex(0) } }, 1, 2), 3, .{}));
    try std.testing.expectEqual(@as(f32, 1000), try penalty(testDof(.{ .scale_fragment = .{} }, 3, 5), 0, .{}));
    try std.testing.expectEqual(@as(f32, 100), try penalty(testDof(.{ .invert_bond = .{ .pivot = core.ids.AtomId.fromIndex(0), .bound = core.ids.AtomId.fromIndex(1) } }, 1, 2), 1, .{}));
    try std.testing.expectEqual(@as(f32, 600), try penalty(testDof(.{ .flip_ring = .{
        .ring = core.ids.RingId.fromIndex(0),
        .pivot_a = core.ids.AtomId.fromIndex(0),
        .pivot_b = core.ids.AtomId.fromIndex(1),
        .penalty_multiplier = 3,
    } }, 1, 2), 0, .{}));
}
