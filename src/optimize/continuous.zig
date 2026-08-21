const std = @import("std");
const core = @import("core");
const geometry = @import("geometry");
const model = @import("model");

pub const marching_squares = @import("marching_squares.zig");
pub const discrete = @import("discrete.zig");

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

pub const ConstructionOptions = struct {
    intrafragment_clashes: bool = true,
    rigid_atoms: []const bool = &.{},
    atom_fragments: []const u32 = &.{},
    atom_has_dofs: []const bool = &.{},
    bond_in_non_macrocycle_ring: []const bool = &.{},
};

pub const BendRingContext = union(enum) {
    non_ring,
    small_ring: u32,
    macrocycle_ring,
};

pub const BendCandidate = struct {
    atom_a: core.ids.AtomId,
    atom_b: core.ids.AtomId,
    initial_rest_degrees: f32 = 120,
    ring: BendRingContext = .non_ring,
};

pub const BendGroup = struct {
    center: core.ids.AtomId,
    /// Adjacent neighbor pairs in clockwise order.
    candidates: []const BendCandidate,
    cross_layout: bool = false,
    inverted_macrocycle_bond: bool = false,
};

pub const BendConstructionOptions = struct {
    even_angles: bool = false,
    rigid_atoms: []const bool = &.{},
};

/// Construct the clash-then-stretch prefix of addInteractionsOfMolecule in
/// upstream insertion order. Bend and E/Z construction depend on ordered ring
/// context and are appended by their owning layout/macrocycle phases.
pub fn buildBaseInteractions(
    allocator: std.mem.Allocator,
    atoms: []const model.Atom,
    bonds: []const model.Bond,
    options: ConstructionOptions,
) core.errors.Error!core.interaction.Collection {
    if ((options.rigid_atoms.len != 0 and options.rigid_atoms.len != atoms.len) or
        (options.atom_fragments.len != 0 and options.atom_fragments.len != atoms.len) or
        (options.atom_has_dofs.len != 0 and options.atom_has_dofs.len != atoms.len) or
        (options.bond_in_non_macrocycle_ring.len != 0 and options.bond_in_non_macrocycle_ring.len != bonds.len))
    {
        return error.InvalidMapping;
    }
    var interactions: std.ArrayList(core.interaction.Interaction) = .empty;
    defer interactions.deinit(allocator);

    if (atoms.len > 1) {
        for (atoms, 0..) |point, point_index| {
            for (bonds) |bond| {
                try validateBond(bond, atoms.len);
                if (point.id == bond.start or point.id == bond.end) continue;
                if (point.fixed and atoms[bond.start.index()].fixed and atoms[bond.end.index()].fixed) continue;
                if (isRigid(options, point_index) and isRigid(options, bond.start.index()) and isRigid(options, bond.end.index())) continue;
                if (areNeighbors(point.id, bond.start, bonds) or areNeighbors(point.id, bond.end, bonds)) continue;
                if (!options.intrafragment_clashes and skipFragmentClash(options, point_index, bond)) continue;

                var rest_scale: f32 = 0.8;
                if (point.atomic_number == .carbon and point.formal_charge == 0) rest_scale -= 0.1;
                const start = atoms[bond.start.index()];
                const end = atoms[bond.end.index()];
                if (start.atomic_number == .carbon and start.formal_charge == 0 and
                    end.atomic_number == .carbon and end.formal_charge == 0)
                {
                    rest_scale -= 0.1;
                }
                const rest = core.math.bond_length * rest_scale;
                try appendInteraction(allocator, &interactions, .{ .clash = .{
                    .segment_start = bond.start,
                    .point = point.id,
                    .segment_end = bond.end,
                    .rest_squared_distance = rest * rest,
                } });
            }
        }
    }

    for (bonds, 0..) |bond, bond_index| {
        try validateBond(bond, atoms.len);
        var force_constant: f32 = 0.1;
        var rest_length = core.math.bond_length;
        if (isRigid(options, bond.start.index()) and isRigid(options, bond.end.index())) {
            rest_length = geometry.length(geometry.subtract(atoms[bond.end.index()].coordinates, atoms[bond.start.index()].coordinates));
        }
        if (options.bond_in_non_macrocycle_ring.len != 0 and options.bond_in_non_macrocycle_ring[bond_index]) {
            force_constant *= 50;
        }
        try appendInteraction(allocator, &interactions, .{ .stretch = .{
            .atom_a = bond.start,
            .atom_b = bond.end,
            .force_constant = force_constant,
            .rest_length = rest_length,
        } });
    }

    return .{
        .allocator = allocator,
        .items = interactions.toOwnedSlice(allocator) catch return error.OutOfMemory,
    };
}

/// Construct bends after topology has supplied clockwise adjacent-neighbor
/// pairs and ring classification. This keeps optimize independent of topology
/// while retaining upstream's group-wide target redistribution and ordering.
pub fn buildBendInteractions(
    allocator: std.mem.Allocator,
    atoms: []const model.Atom,
    groups: []const BendGroup,
    options: BendConstructionOptions,
) core.errors.Error!core.interaction.Collection {
    if (options.rigid_atoms.len != 0 and options.rigid_atoms.len != atoms.len) return error.InvalidMapping;
    var interactions: std.ArrayList(core.interaction.Interaction) = .empty;
    defer interactions.deinit(allocator);
    var bends: std.ArrayList(core.interaction.Bend) = .empty;
    defer bends.deinit(allocator);
    var ring_flags: std.ArrayList(bool) = .empty;
    defer ring_flags.deinit(allocator);

    for (groups) |group| {
        if (!group.center.isValid() or group.center.index() >= atoms.len) return error.InvalidAtomIndex;
        bends.clearRetainingCapacity();
        ring_flags.clearRetainingCapacity();
        var ring_count: usize = 0;
        var non_ring_count: usize = 0;
        var ring_total: f32 = 0;

        for (group.candidates) |candidate| {
            if (!candidate.atom_a.isValid() or candidate.atom_a.index() >= atoms.len or
                !candidate.atom_b.isValid() or candidate.atom_b.index() >= atoms.len or
                candidate.atom_a == group.center or candidate.atom_b == group.center or candidate.atom_a == candidate.atom_b)
            {
                return error.InvalidAtomIndex;
            }
            if (!std.math.isFinite(candidate.initial_rest_degrees)) return error.InvalidOption;
            var bend = core.interaction.Bend{
                .atom_a = candidate.atom_a,
                .center = group.center,
                .atom_b = candidate.atom_b,
                .rest_degrees = candidate.initial_rest_degrees,
            };
            const ring_target = switch (candidate.ring) {
                .non_ring => false,
                .macrocycle_ring => true,
                .small_ring => |effective_size| target: {
                    if (effective_size == 0) return error.InvalidMapping;
                    bend.is_ring = true;
                    bend.force_constant *= 100;
                    bend.rest_degrees = 180 - 360 / @as(f32, @floatFromInt(effective_size));
                    break :target true;
                },
            };
            if (isRigidBend(options, candidate.atom_a.index()) and
                isRigidBend(options, group.center.index()) and
                isRigidBend(options, candidate.atom_b.index()))
            {
                bend.rest_degrees = geometry.unsignedAngle(
                    atoms[candidate.atom_a.index()].coordinates,
                    atoms[group.center.index()].coordinates,
                    atoms[candidate.atom_b.index()].coordinates,
                );
            }
            bends.append(allocator, bend) catch return error.OutOfMemory;
            ring_flags.append(allocator, ring_target) catch return error.OutOfMemory;
            if (ring_target) {
                ring_count += 1;
                ring_total += bend.rest_degrees;
            } else {
                non_ring_count += 1;
            }
        }

        const inverted = group.inverted_macrocycle_bond and ring_count == 1 and non_ring_count == 2;
        if (ring_count != 0) {
            if (inverted) ring_total = 360 - ring_total;
            for (bends.items, ring_flags.items) |*bend, is_ring_target| {
                if (!is_ring_target) bend.rest_degrees = (360 - ring_total) / @as(f32, @floatFromInt(non_ring_count));
            }
        } else if (non_ring_count == 4) {
            if (group.cross_layout or options.even_angles) {
                for (bends.items) |*bend| bend.rest_degrees = 90;
            } else {
                var biggest_index: usize = 0;
                var biggest_angle: f32 = 0;
                for (bends.items, 0..) |bend, index| {
                    const angle = geometry.unsignedAngle(
                        atoms[bend.atom_a.index()].coordinates,
                        atoms[bend.center.index()].coordinates,
                        atoms[bend.atom_b.index()].coordinates,
                    );
                    if (angle > biggest_angle) {
                        biggest_angle = angle;
                        biggest_index = index;
                    }
                }
                bends.items[biggest_index].rest_degrees = 120;
                bends.items[(biggest_index + 1) % 4].rest_degrees = 90;
                bends.items[(biggest_index + 2) % 4].rest_degrees = 60;
                bends.items[(biggest_index + 3) % 4].rest_degrees = 90;
            }
        } else if (non_ring_count > 4) {
            const target = 360 / @as(f32, @floatFromInt(non_ring_count));
            for (bends.items) |*bend| bend.rest_degrees = target;
        }

        for (bends.items) |bend| {
            if (atoms[bend.atom_a.index()].fixed and atoms[bend.center.index()].fixed and atoms[bend.atom_b.index()].fixed) continue;
            try appendInteraction(allocator, &interactions, .{ .bend = bend });
        }
    }

    return .{
        .allocator = allocator,
        .items = interactions.toOwnedSlice(allocator) catch return error.OutOfMemory,
    };
}

pub fn buildConstraintInteractions(
    allocator: std.mem.Allocator,
    atoms: []const model.Atom,
) core.errors.Error!core.interaction.Collection {
    var interactions: std.ArrayList(core.interaction.Interaction) = .empty;
    defer interactions.deinit(allocator);
    for (atoms) |atom| {
        if (!atom.constrained) continue;
        try appendInteraction(allocator, &interactions, .{ .constraint = .{
            .atom = atom.id,
            .origin = atom.template_coordinates orelse return error.InvalidMapping,
        } });
    }
    return .{
        .allocator = allocator,
        .items = interactions.toOwnedSlice(allocator) catch return error.OutOfMemory,
    };
}

/// Materialize E/Z constraints prepared by peptide or macrocycle topology.
pub fn buildEzInteractions(
    allocator: std.mem.Allocator,
    constraints: []const core.interaction.EzConstraint,
) core.errors.Error!core.interaction.Collection {
    var interactions: std.ArrayList(core.interaction.Interaction) = .empty;
    defer interactions.deinit(allocator);
    for (constraints) |constraint| try appendInteraction(allocator, &interactions, .{ .ez_constraint = constraint });
    return .{
        .allocator = allocator,
        .items = interactions.toOwnedSlice(allocator) catch return error.OutOfMemory,
    };
}

/// Join independently constructed phases and assign final insertion-order IDs.
pub fn combineInteractions(
    allocator: std.mem.Allocator,
    collections: []const []const core.interaction.Interaction,
) core.errors.Error!core.interaction.Collection {
    var count: usize = 0;
    for (collections) |collection| {
        count = std.math.add(usize, count, collection.len) catch return error.TooManyItems;
    }
    if (count > std.math.maxInt(u32)) return error.TooManyItems;
    const items = allocator.alloc(core.interaction.Interaction, count) catch return error.OutOfMemory;
    var index: usize = 0;
    for (collections) |collection| {
        for (collection) |interaction| {
            items[index] = interaction;
            items[index].id = core.ids.InteractionId.fromIndex(@intCast(index));
            index += 1;
        }
    }
    return .{ .allocator = allocator, .items = items };
}

fn appendInteraction(
    allocator: std.mem.Allocator,
    interactions: *std.ArrayList(core.interaction.Interaction),
    payload: core.interaction.Payload,
) core.errors.Error!void {
    if (interactions.items.len >= std.math.maxInt(u32)) return error.TooManyItems;
    interactions.append(allocator, .{
        .id = core.ids.InteractionId.fromIndex(@intCast(interactions.items.len)),
        .payload = payload,
    }) catch return error.OutOfMemory;
}

fn validateBond(bond: model.Bond, atom_count: usize) core.errors.Error!void {
    if (!bond.start.isValid() or !bond.end.isValid() or
        bond.start.index() >= atom_count or bond.end.index() >= atom_count or bond.start == bond.end)
    {
        return error.InvalidAtomIndex;
    }
}

fn isRigid(options: ConstructionOptions, index: usize) bool {
    return options.rigid_atoms.len != 0 and options.rigid_atoms[index];
}

fn isRigidBend(options: BendConstructionOptions, index: usize) bool {
    return options.rigid_atoms.len != 0 and options.rigid_atoms[index];
}

fn areNeighbors(first: core.ids.AtomId, second: core.ids.AtomId, bonds: []const model.Bond) bool {
    for (bonds) |bond| {
        if ((bond.start == first and bond.end == second) or (bond.start == second and bond.end == first)) return true;
    }
    return false;
}

fn skipFragmentClash(options: ConstructionOptions, point_index: usize, bond: model.Bond) bool {
    if (options.atom_fragments.len == 0 or options.atom_has_dofs.len == 0) return false;
    const start = bond.start.index();
    const end = bond.end.index();
    if (options.atom_has_dofs[start] or options.atom_has_dofs[point_index] or options.atom_has_dofs[end]) return false;
    return options.atom_fragments[start] == options.atom_fragments[point_index] or
        options.atom_fragments[end] == options.atom_fragments[point_index];
}

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

fn baseInteractionsAndDiscard(allocator: std.mem.Allocator) !void {
    const atoms = [_]model.Atom{
        .{ .id = core.ids.AtomId.fromIndex(0), .input_index = 0, .atomic_number = .carbon },
        .{ .id = core.ids.AtomId.fromIndex(1), .input_index = 1, .atomic_number = .carbon, .coordinates = .{ .x = 60 } },
        .{ .id = core.ids.AtomId.fromIndex(2), .input_index = 2, .atomic_number = .oxygen, .coordinates = .{ .x = 30, .y = 5 } },
    };
    const bonds = [_]model.Bond{.{
        .id = core.ids.BondId.fromIndex(0),
        .input_index = 0,
        .start = atoms[0].id,
        .end = atoms[1].id,
        .input_order = .single,
        .effective_order = .single,
    }};
    var interactions = try buildBaseInteractions(allocator, &atoms, &bonds, .{
        .rigid_atoms = &.{ true, true, false },
        .bond_in_non_macrocycle_ring = &.{true},
    });
    defer interactions.deinit();
    try std.testing.expectEqual(@as(usize, 2), interactions.items.len);
    const clash = interactions.items[0].payload.clash;
    try std.testing.expectEqual(@as(f32, 1225), clash.rest_squared_distance);
    const stretch = interactions.items[1].payload.stretch;
    try std.testing.expectEqual(@as(f32, 5), stretch.force_constant);
    try std.testing.expectEqual(@as(f32, 60), stretch.rest_length);
}

test "base interaction construction preserves clash-stretch order and cleans allocation failures" {
    try std.testing.checkAllAllocationFailures(std.testing.allocator, baseInteractionsAndDiscard, .{});
}

fn bendInteractionsAndDiscard(allocator: std.mem.Allocator) !void {
    const atoms = [_]model.Atom{
        .{ .id = core.ids.AtomId.fromIndex(0), .input_index = 0, .atomic_number = .carbon },
        .{ .id = core.ids.AtomId.fromIndex(1), .input_index = 1, .atomic_number = .carbon, .coordinates = .{ .x = 1 } },
        .{ .id = core.ids.AtomId.fromIndex(2), .input_index = 2, .atomic_number = .carbon, .coordinates = .{ .y = 1 } },
        .{ .id = core.ids.AtomId.fromIndex(3), .input_index = 3, .atomic_number = .carbon, .coordinates = .{ .x = -1 } },
    };
    const candidates = [_]BendCandidate{
        .{ .atom_a = atoms[1].id, .atom_b = atoms[2].id, .ring = .{ .small_ring = 4 } },
        .{ .atom_a = atoms[2].id, .atom_b = atoms[3].id },
        .{ .atom_a = atoms[3].id, .atom_b = atoms[1].id },
    };
    var interactions = try buildBendInteractions(allocator, &atoms, &.{.{
        .center = atoms[0].id,
        .candidates = &candidates,
    }}, .{});
    defer interactions.deinit();
    try std.testing.expectEqual(@as(usize, 3), interactions.items.len);
    try std.testing.expectEqual(@as(f32, 100), interactions.items[0].payload.bend.force_constant);
    try std.testing.expectEqual(@as(f32, 90), interactions.items[0].payload.bend.rest_degrees);
    try std.testing.expectEqual(@as(f32, 135), interactions.items[1].payload.bend.rest_degrees);
    try std.testing.expectEqual(@as(f32, 135), interactions.items[2].payload.bend.rest_degrees);
}

test "bend construction redistributes ring angles and cleans allocation failures" {
    try std.testing.checkAllAllocationFailures(std.testing.allocator, bendInteractionsAndDiscard, .{});
}

fn remainingInteractionsAndDiscard(allocator: std.mem.Allocator) !void {
    const atoms = [_]model.Atom{.{
        .id = core.ids.AtomId.fromIndex(0),
        .input_index = 0,
        .atomic_number = .carbon,
        .constrained = true,
        .template_coordinates = .{ .x = 2, .y = 3 },
    }};
    var constraints = try buildConstraintInteractions(allocator, &atoms);
    defer constraints.deinit();
    var ez = try buildEzInteractions(allocator, &.{.{
        .side_a = core.ids.AtomId.fromIndex(0),
        .double_a = core.ids.AtomId.fromIndex(1),
        .double_b = core.ids.AtomId.fromIndex(2),
        .side_b = core.ids.AtomId.fromIndex(3),
        .is_z = true,
    }});
    defer ez.deinit();
    var combined = try combineInteractions(allocator, &.{ constraints.items, ez.items });
    defer combined.deinit();
    try std.testing.expectEqual(@as(usize, 2), combined.items.len);
    try std.testing.expectEqual(core.ids.InteractionId.fromIndex(1), combined.items[1].id);
}

test "constraint and E/Z construction combine with complete allocation cleanup" {
    try std.testing.checkAllAllocationFailures(std.testing.allocator, remainingInteractionsAndDiscard, .{});
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
    // Same-process pinned C++ probe: 4.99363375, 55.0064354.
    try std.testing.expectApproxEqAbs(@as(f32, 4.99363375), coordinates[0].x, 0.0002);
    try std.testing.expectApproxEqAbs(@as(f32, 55.0064354), coordinates[1].x, 0.0002);
    try std.testing.expectApproxEqAbs(@as(f32, 50), geometry.length(geometry.subtract(coordinates[0], coordinates[1])), 0.05);

    const before = coordinates;
    const skipped = try minimize(&interactions, state, .{ .skip = true });
    try std.testing.expectEqual(@as(usize, 0), skipped.iterations);
    try std.testing.expectEqual(before, coordinates);
}

test "fixed-center bend minimization agrees with pinned oracle fixture" {
    var coordinates = [_]core.math.Vec2{ .{ .x = 50 }, .{}, .{ .y = 50 } };
    var forces = [_]core.math.Vec2{ .{}, .{}, .{} };
    const fixed = [_]bool{ false, true, false };
    const interactions = [_]core.interaction.Interaction{.{
        .id = core.ids.InteractionId.fromIndex(0),
        .payload = .{ .bend = .{
            .atom_a = core.ids.AtomId.fromIndex(0),
            .center = core.ids.AtomId.fromIndex(1),
            .atom_b = core.ids.AtomId.fromIndex(2),
        } },
    }};
    _ = try minimize(&interactions, .{ .coordinates = &coordinates, .forces = &forces, .fixed = &fixed }, .{});
    // Same-process pinned C++ probe with the center constrained/fixed.
    const expected = [_]core.math.Vec2{
        .{ .x = 48.4643898, .y = -12.3434944 },
        .{},
        .{ .x = -12.3434944, .y = 48.4643898 },
    };
    for (expected, coordinates) |want, actual| {
        try std.testing.expectApproxEqAbs(want.x, actual.x, 0.001);
        try std.testing.expectApproxEqAbs(want.y, actual.y, 0.001);
    }
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
