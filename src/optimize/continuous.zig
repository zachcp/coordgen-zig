const std = @import("std");
const core = @import("core");
const geometry = @import("geometry");
const model = @import("model");

pub const marching_squares = @import("marching_squares.zig");

pub const State = struct {
    coordinates: []core.math.Vec2,
    forces: []core.math.Vec2,
    fixed: []const bool,

    fn validate(self: State) core.errors.Error!void {
        if (self.coordinates.len != self.forces.len or self.coordinates.len != self.fixed.len) return error.InvalidMapping;
    }
};

pub fn score(interaction: core.interaction.Interaction, state: State) core.errors.Error!f32 {
    try state.validate();
    return switch (interaction.payload) {
        .stretch => |value| scoreStretch(value, state),
        .bend => |value| scoreBend(value, state),
        .clash => |value| scoreClash(value, state, false),
        .constraint => |value| scoreConstraint(value, state),
        .ez_constraint => |value| scoreEz(value, state),
    };
}

pub fn scoreAll(interactions: []const core.interaction.Interaction, state: State) core.errors.Error!f32 {
    try state.validate();
    var energy: f32 = 0;
    for (interactions) |interaction| energy += try score(interaction, state);
    return energy;
}

pub const MinimizeOptions = struct {
    skip: bool = false,
    max_iterations: usize = 1000,
    max_step: f32 = 0.1,
};

pub const MinimizeResult = struct {
    iterations: usize,
    energy: f32,
    converged: bool,
};

pub const StereoValidator = *const fn ([]const model.Atom) bool;

/// Run the pinned continuous minimization loop. Interaction scoring both
/// returns energy and accumulates force; movable force is reset by applyForces.
pub fn minimize(
    interactions: []const core.interaction.Interaction,
    state: State,
    options: MinimizeOptions,
) core.errors.Error!MinimizeResult {
    try state.validate();
    if (!std.math.isFinite(options.max_step) or options.max_step < 0 or options.max_iterations > 1000) {
        return error.InvalidOption;
    }
    if (options.skip or options.max_iterations == 0) return .{ .iterations = 0, .energy = 0, .converged = true };

    var energies: [1000]f32 = undefined;
    for (0..options.max_iterations) |iteration| {
        const energy = try scoreAll(interactions, state);
        energies[iteration] = energy;
        if (!try applyForces(state, options.max_step)) {
            return .{ .iterations = iteration + 1, .energy = energy, .converged = true };
        }
        if (iteration >= 200 and energies[iteration - 100] - energy < 20) {
            return .{ .iterations = iteration + 1, .energy = energy, .converged = false };
        }
    }
    return .{
        .iterations = options.max_iterations,
        .energy = energies[options.max_iterations - 1],
        .converged = false,
    };
}

/// Adapt model atoms to the structure-of-arrays minimizer and preserve the
/// upstream molecule rollback when final stereochemistry validation fails.
pub fn minimizeMolecule(
    allocator: std.mem.Allocator,
    atoms: []model.Atom,
    interactions: []const core.interaction.Interaction,
    options: MinimizeOptions,
    validate_stereo: ?StereoValidator,
) core.errors.Error!MinimizeResult {
    const coordinates = allocator.alloc(core.math.Vec2, atoms.len) catch return error.OutOfMemory;
    defer allocator.free(coordinates);
    const previous_coordinates = allocator.alloc(core.math.Vec2, atoms.len) catch return error.OutOfMemory;
    defer allocator.free(previous_coordinates);
    const forces = allocator.alloc(core.math.Vec2, atoms.len) catch return error.OutOfMemory;
    defer allocator.free(forces);
    const fixed = allocator.alloc(bool, atoms.len) catch return error.OutOfMemory;
    defer allocator.free(fixed);
    for (atoms, coordinates, previous_coordinates, forces, fixed) |atom, *position, *previous, *force, *is_fixed| {
        position.* = atom.coordinates;
        previous.* = atom.coordinates;
        force.* = .{};
        is_fixed.* = atom.fixed;
    }

    const result = try minimize(interactions, .{
        .coordinates = coordinates,
        .forces = forces,
        .fixed = fixed,
    }, options);
    for (atoms, coordinates) |*atom, position| atom.coordinates = position;
    if (validate_stereo) |validator| {
        if (!validator(atoms)) {
            for (atoms, previous_coordinates) |*atom, previous| {
                atom.coordinates = previous;
            }
        }
    }
    return result;
}

pub fn hasNaNCoordinates(coordinates: []const core.math.Vec2) bool {
    for (coordinates) |position| {
        if (std.math.isNan(position.x) or std.math.isNan(position.y)) return true;
    }
    return false;
}

pub fn hasValid3DCoordinates(coordinates: []const ?core.math.Vec3) bool {
    const invalid_coordinates: f32 = 10_000_001;
    for (coordinates) |optional| {
        const position = optional orelse return false;
        if (!(position.x < invalid_coordinates and position.y < invalid_coordinates and position.z < invalid_coordinates)) return false;
    }
    return true;
}

/// Apply the upstream emergency x/-y projection and two-decimal rounding.
pub fn fallbackOn3DCoordinates(
    coordinates: []core.math.Vec2,
    coordinates_3d: []const ?core.math.Vec3,
) core.errors.Error!void {
    if (coordinates.len != coordinates_3d.len or !hasValid3DCoordinates(coordinates_3d)) return error.InvalidCoordinate;
    for (coordinates, coordinates_3d) |*position, optional| {
        const source = optional.?;
        position.* = .{
            .x = roundToTwoDecimalDigits(source.x * 35),
            .y = roundToTwoDecimalDigits(-source.y * 35),
        };
    }
}

fn roundToTwoDecimalDigits(value: f32) f32 {
    return @floor(value * 100 + 0.5) * 0.01;
}

fn coordinate(state: State, atom: core.ids.AtomId) core.errors.Error!core.math.Vec2 {
    if (!atom.isValid() or atom.index() >= state.coordinates.len) return error.InvalidAtomIndex;
    return state.coordinates[atom.index()];
}

fn addForce(state: State, atom: core.ids.AtomId, force: core.math.Vec2) void {
    state.forces[atom.index()] = geometry.add(state.forces[atom.index()], force);
}

fn scoreStretch(value: core.interaction.Stretch, state: State) core.errors.Error!f32 {
    const difference = geometry.subtract(try coordinate(state, value.atom_a), try coordinate(state, value.atom_b));
    const magnitude = geometry.length(difference);
    const delta = magnitude - value.rest_length;
    const energy = 0.5 * value.force_constant * delta * delta;
    const force_delta = if (magnitude < value.rest_length - value.tolerance)
        value.rest_length - value.tolerance - magnitude
    else if (magnitude > value.rest_length + value.tolerance)
        value.rest_length + value.tolerance - magnitude
    else
        return energy;
    const short_penalty = @max(@as(f32, 0), 0.4 * value.rest_length - magnitude);
    const force = geometry.scale(geometry.normalize(difference), value.force_constant * force_delta + 10 * short_penalty);
    addForce(state, value.atom_a, force);
    addForce(state, value.atom_b, geometry.negate(force));
    return energy;
}

fn scoreBend(value: core.interaction.Bend, state: State) core.errors.Error!f32 {
    const first = try coordinate(state, value.atom_a);
    const center = try coordinate(state, value.center);
    const last = try coordinate(state, value.atom_b);
    const angle = geometry.unsignedAngle(first, center, last);
    const energy_delta = angle - value.rest_degrees;
    const energy = 5 * value.force_constant * value.secondary_force_constant * energy_delta * energy_delta;
    const target = if (value.rest_degrees > 180) 360 - value.rest_degrees else value.rest_degrees;
    const force_delta = target - @abs(angle);
    const arm_a = geometry.subtract(first, center);
    const arm_b = geometry.subtract(last, center);
    const chord = geometry.subtract(last, first);
    var normal_a = core.math.Vec2{ .x = arm_a.y, .y = -arm_a.x };
    var normal_b = core.math.Vec2{ .x = arm_b.y, .y = -arm_b.x };
    if (geometry.dot(normal_a, chord) > 0) normal_a = geometry.negate(normal_a);
    if (geometry.dot(normal_b, chord) < 0) normal_b = geometry.negate(normal_b);
    normal_a = geometry.scale(normalizeWithFloor(normal_a), value.force_constant * value.secondary_force_constant * force_delta);
    normal_b = geometry.scale(normalizeWithFloor(normal_b), value.force_constant * value.secondary_force_constant * force_delta);
    addForce(state, value.atom_a, normal_a);
    addForce(state, value.atom_b, normal_b);
    addForce(state, value.center, geometry.negate(geometry.add(normal_a, normal_b)));
    return energy;
}

fn normalizeWithFloor(value: core.math.Vec2) core.math.Vec2 {
    return geometry.divide(value, @max(geometry.length(value), geometry.epsilon));
}

fn scoreClash(value: core.interaction.Clash, state: State, skip_force: bool) core.errors.Error!f32 {
    const start = try coordinate(state, value.segment_start);
    const point = try coordinate(state, value.point);
    const end = try coordinate(state, value.segment_end);
    const distance = geometry.squaredDistancePointSegment(point, start, end).squared_distance;
    if (distance > value.rest_squared_distance) return 0;
    const deficit = value.rest_squared_distance - distance;
    const energy = 0.5 * value.force_constant * value.secondary_force_constant * deficit;
    if (skip_force) return energy;
    const projection = geometry.projectPointOnLine(point, start, end);
    const force = geometry.scale(geometry.normalize(geometry.subtract(point, projection)), deficit * value.force_constant * value.secondary_force_constant);
    addForce(state, value.point, force);
    addForce(state, value.segment_start, geometry.scale(force, -0.5));
    addForce(state, value.segment_end, geometry.scale(force, -0.5));
    return energy;
}

fn scoreConstraint(value: core.interaction.Constraint, state: State) core.errors.Error!f32 {
    return value.force_constant * geometry.squaredDistance(try coordinate(state, value.atom), value.origin);
}

fn scoreEz(value: core.interaction.EzConstraint, state: State) core.errors.Error!f32 {
    const side_a = try coordinate(state, value.side_a);
    const double_a = try coordinate(state, value.double_a);
    const double_b = try coordinate(state, value.double_b);
    const side_b = try coordinate(state, value.side_b);
    if (geometry.sameSide(side_a, side_b, double_a, double_b) == value.is_z) return 0;
    const projection_a = geometry.projectPointOnLine(side_a, double_a, double_b);
    const projection_b = geometry.projectPointOnLine(side_b, double_a, double_b);
    const choose_b = geometry.squaredDistance(side_a, projection_a) > geometry.squaredDistance(side_b, projection_b);
    const side = if (choose_b) value.side_b else value.side_a;
    const center = if (choose_b) value.double_b else value.double_a;
    var force = geometry.subtract(if (choose_b) projection_b else projection_a, if (choose_b) side_b else side_a);
    if (value.force_movement) {
        state.coordinates[side.index()] = geometry.add(state.coordinates[side.index()], force);
        state.coordinates[center.index()] = geometry.subtract(state.coordinates[center.index()], force);
        state.forces[side.index()] = .{};
        state.forces[center.index()] = .{};
    } else {
        force = geometry.scale(geometry.normalize(force), 10);
        addForce(state, side, force);
        addForce(state, center, geometry.negate(force));
    }
    return 5000;
}

pub fn applyForces(state: State, max_step: f32) core.errors.Error!bool {
    try state.validate();
    if (!std.math.isFinite(max_step) or max_step < 0) return error.InvalidOption;
    var total_squared_displacement: f32 = 0;
    for (state.coordinates, state.forces, state.fixed) |*position, *force, fixed| {
        if (fixed) continue;
        var displacement = geometry.scale(force.*, 0.3);
        if (std.math.isNan(displacement.x) or std.math.isNan(displacement.y)) displacement = .{};
        var squared = geometry.squaredLength(displacement);
        if (squared < geometry.epsilon) squared = geometry.epsilon;
        if (squared > max_step * max_step) displacement = geometry.scale(displacement, max_step / @sqrt(squared));
        position.* = geometry.add(position.*, displacement);
        total_squared_displacement += geometry.squaredLength(displacement);
        force.* = .{};
    }
    return total_squared_displacement >= 0.001;
}

test "all continuous interactions conserve direct energy and force facts" {
    var coordinates = [_]core.math.Vec2{ .{}, .{ .x = 60 }, .{ .x = 60, .y = 40 }, .{ .x = 120 } };
    var forces = [_]core.math.Vec2{ .{}, .{}, .{}, .{} };
    const fixed = [_]bool{ false, false, false, false };
    const state = State{ .coordinates = &coordinates, .forces = &forces, .fixed = &fixed };
    const stretch = core.interaction.Interaction{ .id = core.ids.InteractionId.fromIndex(0), .payload = .{ .stretch = .{
        .atom_a = core.ids.AtomId.fromIndex(0),
        .atom_b = core.ids.AtomId.fromIndex(1),
    } } };
    try std.testing.expectEqual(@as(f32, 50), try score(stretch, state));
    try std.testing.expectApproxEqAbs(@as(f32, 10), forces[0].x, 0.001);
    try std.testing.expectApproxEqAbs(@as(f32, -10), forces[1].x, 0.001);
    @memset(&forces, .{});
    const constraint = core.interaction.Interaction{ .id = core.ids.InteractionId.fromIndex(1), .payload = .{ .constraint = .{
        .atom = core.ids.AtomId.fromIndex(1),
        .origin = .{ .x = 50 },
    } } };
    try std.testing.expectEqual(@as(f32, 50), try score(constraint, state));
    try std.testing.expectEqual(core.math.Vec2{}, forces[1]);
}

test "bend clash and E/Z values and forces match pinned formulas" {
    var coordinates = [_]core.math.Vec2{
        .{ .x = 1 },
        .{},
        .{ .y = 1 },
        .{ .x = -1 },
        .{ .x = 1, .y = 1 },
    };
    var forces = [_]core.math.Vec2{ .{}, .{}, .{}, .{}, .{} };
    const fixed = [_]bool{ false, false, false, false, false };
    const state = State{ .coordinates = &coordinates, .forces = &forces, .fixed = &fixed };

    const bend = core.interaction.Interaction{ .id = core.ids.InteractionId.fromIndex(0), .payload = .{ .bend = .{
        .atom_a = core.ids.AtomId.fromIndex(0),
        .center = core.ids.AtomId.fromIndex(1),
        .atom_b = core.ids.AtomId.fromIndex(2),
    } } };
    try std.testing.expectApproxEqAbs(@as(f32, 225), try score(bend, state), 0.001);
    try std.testing.expectApproxEqAbs(@as(f32, -1.5), forces[0].y, 0.001);
    try std.testing.expectApproxEqAbs(@as(f32, -1.5), forces[2].x, 0.001);
    try expectForceSumZero(&forces);

    @memset(&forces, .{});
    const clash = core.interaction.Interaction{ .id = core.ids.InteractionId.fromIndex(1), .payload = .{ .clash = .{
        .segment_start = core.ids.AtomId.fromIndex(3),
        .point = core.ids.AtomId.fromIndex(2),
        .segment_end = core.ids.AtomId.fromIndex(0),
    } } };
    try std.testing.expectApproxEqAbs(@as(f32, 44.95), try score(clash, state), 0.001);
    try std.testing.expectApproxEqAbs(@as(f32, 89.9), forces[2].y, 0.001);
    try expectForceSumZero(&forces);

    @memset(&forces, .{});
    const ez = core.interaction.Interaction{ .id = core.ids.InteractionId.fromIndex(2), .payload = .{ .ez_constraint = .{
        .side_a = core.ids.AtomId.fromIndex(2),
        .double_a = core.ids.AtomId.fromIndex(1),
        .double_b = core.ids.AtomId.fromIndex(0),
        .side_b = core.ids.AtomId.fromIndex(4),
        .is_z = false,
    } } };
    try std.testing.expectEqual(@as(f32, 5000), try score(ez, state));
    try std.testing.expectEqual(@as(f32, -10), forces[2].y);
    try std.testing.expectEqual(@as(f32, 10), forces[1].y);
    try expectForceSumZero(&forces);
}

fn expectForceSumZero(forces: []const core.math.Vec2) !void {
    var sum: core.math.Vec2 = .{};
    for (forces) |force| sum = geometry.add(sum, force);
    try std.testing.expectApproxEqAbs(@as(f32, 0), sum.x, 0.0001);
    try std.testing.expectApproxEqAbs(@as(f32, 0), sum.y, 0.0001);
}

test "force application caps movement, clears movable force, and retains fixed force" {
    var coordinates = [_]core.math.Vec2{ .{}, .{} };
    var forces = [_]core.math.Vec2{ .{ .x = 10 }, .{ .x = 10 } };
    const fixed = [_]bool{ false, true };
    try std.testing.expect(try applyForces(.{ .coordinates = &coordinates, .forces = &forces, .fixed = &fixed }, 0.1));
    try std.testing.expectApproxEqAbs(@as(f32, 0.1), coordinates[0].x, 0.0001);
    try std.testing.expectEqual(@as(f32, 10), forces[1].x);
}

test "continuous minimization converges a stretched bond and honors skip" {
    var coordinates = [_]core.math.Vec2{ .{}, .{ .x = 60 } };
    var forces = [_]core.math.Vec2{ .{}, .{} };
    const fixed = [_]bool{ false, false };
    const interactions = [_]core.interaction.Interaction{.{
        .id = core.ids.InteractionId.fromIndex(0),
        .payload = .{ .stretch = .{
            .atom_a = core.ids.AtomId.fromIndex(0),
            .atom_b = core.ids.AtomId.fromIndex(1),
        } },
    }};
    const state = State{ .coordinates = &coordinates, .forces = &forces, .fixed = &fixed };
    const result = try minimize(&interactions, state, .{});
    try std.testing.expect(result.converged);
    try std.testing.expect(result.iterations > 1);
    try std.testing.expectApproxEqAbs(@as(f32, 50), geometry.length(geometry.subtract(coordinates[0], coordinates[1])), 0.05);

    const before = coordinates;
    const skipped = try minimize(&interactions, state, .{ .skip = true });
    try std.testing.expectEqual(@as(usize, 0), skipped.iterations);
    try std.testing.expectEqual(before, coordinates);
}

test "nonfinite detection and 3D fallback preserve pinned asymmetries" {
    try std.testing.expect(hasNaNCoordinates(&.{.{ .x = std.math.nan(f32) }}));
    try std.testing.expect(!hasNaNCoordinates(&.{.{ .x = std.math.inf(f32) }}));
    const source = [_]?core.math.Vec3{.{ .x = 1.234, .y = -2.345, .z = 7 }};
    try std.testing.expect(hasValid3DCoordinates(&source));
    var coordinates = [_]core.math.Vec2{.{}};
    try fallbackOn3DCoordinates(&coordinates, &source);
    try std.testing.expectEqual(core.math.Vec2{ .x = 43.19, .y = 82.08 }, coordinates[0]);
    try std.testing.expect(hasValid3DCoordinates(&.{.{ .x = -std.math.inf(f32) }}));
    try std.testing.expect(!hasValid3DCoordinates(&.{null}));
}

fn rejectStereo(_: []const model.Atom) bool {
    return false;
}

fn minimizeMoleculeAndDiscard(allocator: std.mem.Allocator) !void {
    var atoms = [_]model.Atom{
        .{ .id = core.ids.AtomId.fromIndex(0), .input_index = 0, .atomic_number = .carbon },
        .{ .id = core.ids.AtomId.fromIndex(1), .input_index = 1, .atomic_number = .carbon, .coordinates = .{ .x = 60 } },
    };
    const interactions = [_]core.interaction.Interaction{.{
        .id = core.ids.InteractionId.fromIndex(0),
        .payload = .{ .stretch = .{ .atom_a = atoms[0].id, .atom_b = atoms[1].id } },
    }};
    _ = try minimizeMolecule(allocator, &atoms, &interactions, .{}, rejectStereo);
    try std.testing.expectEqual(@as(f32, 0), atoms[0].coordinates.x);
    try std.testing.expectEqual(@as(f32, 60), atoms[1].coordinates.x);
}

test "molecule minimization rolls back invalid stereo and cleans every allocation failure" {
    try std.testing.checkAllAllocationFailures(std.testing.allocator, minimizeMoleculeAndDiscard, .{});
}

test {
    _ = marching_squares;
}
