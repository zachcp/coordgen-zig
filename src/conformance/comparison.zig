//! Typed comparison semantics for the future native-vs-oracle runner.
//!
//! A comparison mode is selected for every input *and* observable. Selecting
//! one mode for a whole molecule would hide stable structural evidence merely
//! because that molecule's coordinates are unstable.

const std = @import("std");

pub const Tier = enum { exact, tolerant, invariant_statistical };

pub const Observable = enum {
    status,
    clean_pose,
    coordinates,
    input_to_internal,
    internal_to_input,
    effective_bond_orders,
    bond_displays,
    atom_stereo,
    probe_status,
    probe_clean_pose,
    morgan_ranks,
    rings,
    rings_set,
    fragments,
    fragments_set,
    dofs,
    dofs_set,
    dof_penalties,
    dof_penalties_set,
    template_mappings,
    template_mappings_set,
    components,
    components_set,
    component_transforms,
    component_transforms_set,
};

pub const Stability = struct {
    architecture: bool,
    allocator_order: bool,

    pub fn exact(self: Stability) bool {
        return self.architecture and self.allocator_order;
    }
};

/// Structural values are exact whenever their own stability record is exact.
/// Coordinates and derived floats use normalized/tolerant comparisons only if
/// the oracle is stable, otherwise they use the aggregate quality tier.
pub fn tierFor(observable: Observable, stability: Stability) Tier {
    if (!stability.exact()) return .invariant_statistical;
    return switch (observable) {
        .coordinates,
        .dof_penalties,
        .dof_penalties_set,
        .component_transforms,
        .component_transforms_set,
        => .tolerant,
        else => .exact,
    };
}

pub const bond_length: f64 = 50;

pub const Point = struct { x: f64, y: f64 };

/// Reflection is unavailable unless the fixture records an achirality proof.
/// The default forbids it, preventing a mirrored chiral layout from passing a
/// translation/rotation-normalized comparison.
pub const ReflectionPolicy = struct {
    achiral_proof: ?[]const u8 = null,

    pub fn allowsReflection(self: ReflectionPolicy) bool {
        return if (self.achiral_proof) |proof| proof.len != 0 else false;
    }
};

pub const CoordinateResult = struct {
    matches: bool,
    max_deviation_bond_lengths: f64,
};

/// Compare indexed 2D coordinates after translation and rotation. A second,
/// reflected alignment is attempted only when the fixture opts in with proof.
pub fn compareCoordinates(
    reference: []const Point,
    candidate: []const Point,
    tolerance_bond_lengths: f64,
    reflection: ReflectionPolicy,
) !CoordinateResult {
    if (reference.len != candidate.len) return error.CountMismatch;
    if (!std.math.isFinite(tolerance_bond_lengths) or tolerance_bond_lengths < 0) {
        return error.InvalidTolerance;
    }
    const direct = try alignedDeviation(reference, candidate, false);
    if (direct <= tolerance_bond_lengths) return .{
        .matches = true,
        .max_deviation_bond_lengths = direct,
    };
    if (!reflection.allowsReflection()) return .{
        .matches = false,
        .max_deviation_bond_lengths = direct,
    };
    const reflected_deviation = try alignedDeviation(reference, candidate, true);
    return .{
        .matches = reflected_deviation <= tolerance_bond_lengths,
        .max_deviation_bond_lengths = reflected_deviation,
    };
}

fn alignedDeviation(reference: []const Point, candidate: []const Point, reflect: bool) !f64 {
    if (reference.len == 0) return 0;
    if (reference.len == 1) return distance(reference[0], candidate[0]) / bond_length;

    const candidate_first = reflected(candidate[0], reflect);
    var anchor: ?usize = null;
    for (reference[1..], candidate[1..], 1..) |reference_point, candidate_point, index| {
        if (norm(subtract(reference_point, reference[0])) != 0 and
            norm(subtract(reflected(candidate_point, reflect), candidate_first)) != 0)
        {
            anchor = index;
            break;
        }
    }
    const anchor_index = anchor orelse return error.DegenerateAlignment;
    const reference_vector = subtract(reference[anchor_index], reference[0]);
    const candidate_vector = subtract(reflected(candidate[anchor_index], reflect), candidate_first);
    const denominator = norm(reference_vector) * norm(candidate_vector);
    const cosine = dot(candidate_vector, reference_vector) / denominator;
    const sine = cross(candidate_vector, reference_vector) / denominator;

    var max_deviation: f64 = 0;
    for (reference, candidate) |reference_point, candidate_point| {
        const translated = subtract(reflected(candidate_point, reflect), candidate_first);
        const rotated = Point{
            .x = translated.x * cosine - translated.y * sine + reference[0].x,
            .y = translated.x * sine + translated.y * cosine + reference[0].y,
        };
        max_deviation = @max(max_deviation, distance(reference_point, rotated) / bond_length);
    }
    return max_deviation;
}

fn reflected(point: Point, enabled: bool) Point {
    return if (enabled) .{ .x = point.x, .y = -point.y } else point;
}

fn subtract(left: Point, right: Point) Point {
    return .{ .x = left.x - right.x, .y = left.y - right.y };
}

fn norm(point: Point) f64 {
    return std.math.sqrt(dot(point, point));
}

fn dot(left: Point, right: Point) f64 {
    return left.x * right.x + left.y * right.y;
}

fn cross(left: Point, right: Point) f64 {
    return left.x * right.y - left.y * right.x;
}

fn distance(left: Point, right: Point) f64 {
    return norm(subtract(left, right));
}

/// Coordinate parity is valid only within one process/build/target. This
/// identity belongs to each baseline instead of committing cross-platform
/// coordinate goldens that would create false failures.
pub const RunIdentity = struct {
    process_id: []const u8,
    target: []const u8,
    toolchain: []const u8,
    optimize_mode: []const u8,

    pub fn requireSameBuild(oracle: RunIdentity, native: RunIdentity) !void {
        if (!std.mem.eql(u8, oracle.process_id, native.process_id) or
            !std.mem.eql(u8, oracle.target, native.target) or
            !std.mem.eql(u8, oracle.toolchain, native.toolchain) or
            !std.mem.eql(u8, oracle.optimize_mode, native.optimize_mode))
        {
            return error.MismatchedRunIdentity;
        }
    }
};

/// Lower values are better. The invariant tier aggregates these per corpus
/// partition and requires native to be no worse than oracle plus a named
/// margin for every metric.
pub const QualityMetrics = struct {
    clash_score: f64,
    bond_length_rms: f64,
    bond_angle_deviation: f64,
    bond_crossings: u32,
    atoms_inside_rings: u32,
};

pub const QualityMargins = struct {
    clash_score: f64 = 0,
    bond_length_rms: f64 = 0,
    bond_angle_deviation: f64 = 0,
    bond_crossings: u32 = 0,
    atoms_inside_rings: u32 = 0,
};

pub fn noWorseThan(oracle: QualityMetrics, native: QualityMetrics, margin: QualityMargins) bool {
    return std.math.isFinite(oracle.clash_score) and
        std.math.isFinite(oracle.bond_length_rms) and
        std.math.isFinite(oracle.bond_angle_deviation) and
        std.math.isFinite(native.clash_score) and
        std.math.isFinite(native.bond_length_rms) and
        std.math.isFinite(native.bond_angle_deviation) and
        std.math.isFinite(margin.clash_score) and
        std.math.isFinite(margin.bond_length_rms) and
        std.math.isFinite(margin.bond_angle_deviation) and
        margin.clash_score >= 0 and
        margin.bond_length_rms >= 0 and
        margin.bond_angle_deviation >= 0 and
        native.clash_score <= oracle.clash_score + margin.clash_score and
        native.bond_length_rms <= oracle.bond_length_rms + margin.bond_length_rms and
        native.bond_angle_deviation <= oracle.bond_angle_deviation + margin.bond_angle_deviation and
        native.bond_crossings <= oracle.bond_crossings +| margin.bond_crossings and
        native.atoms_inside_rings <= oracle.atoms_inside_rings +| margin.atoms_inside_rings;
}

test "tier selection stays per input and observable" {
    const stable = Stability{ .architecture = true, .allocator_order = true };
    const unstable_coordinates = Stability{ .architecture = false, .allocator_order = true };
    try std.testing.expectEqual(Tier.exact, tierFor(.rings, stable));
    try std.testing.expectEqual(Tier.tolerant, tierFor(.coordinates, stable));
    try std.testing.expectEqual(Tier.invariant_statistical, tierFor(.coordinates, unstable_coordinates));
}

test "translation and rotation pass while reflection fails by default" {
    const reference = [_]Point{ .{ .x = 0, .y = 0 }, .{ .x = 50, .y = 0 }, .{ .x = 12.5, .y = 37.5 } };
    const rotated = [_]Point{ .{ .x = 200, .y = 300 }, .{ .x = 200, .y = 350 }, .{ .x = 162.5, .y = 312.5 } };
    const mirrored = [_]Point{ .{ .x = 0, .y = 0 }, .{ .x = 50, .y = 0 }, .{ .x = 12.5, .y = -37.5 } };

    try std.testing.expect((try compareCoordinates(&reference, &rotated, 1e-12, .{})).matches);
    try std.testing.expect(!(try compareCoordinates(&reference, &mirrored, 1e-12, .{})).matches);
    try std.testing.expect((try compareCoordinates(&reference, &mirrored, 1e-12, .{
        .achiral_proof = "fixture demonstrated achirality",
    })).matches);
}

test "coordinate comparisons require an identical run identity" {
    const native = RunIdentity{
        .process_id = "run-4",
        .target = "aarch64-macos",
        .toolchain = "0.17.0-dev.1516+8a4b5424d",
        .optimize_mode = "ReleaseFast",
    };
    try RunIdentity.requireSameBuild(native, native);
    var other_target = native;
    other_target.target = "x86_64-macos";
    try std.testing.expectError(error.MismatchedRunIdentity, RunIdentity.requireSameBuild(native, other_target));
}

test "invariant quality floor rejects a worse native result" {
    const oracle = QualityMetrics{
        .clash_score = 1,
        .bond_length_rms = 0.1,
        .bond_angle_deviation = 0.2,
        .bond_crossings = 2,
        .atoms_inside_rings = 0,
    };
    try std.testing.expect(noWorseThan(oracle, oracle, .{}));
    var worse = oracle;
    worse.bond_crossings += 1;
    try std.testing.expect(!noWorseThan(oracle, worse, .{}));
}
