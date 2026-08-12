const std = @import("std");
const core = @import("core");

pub const Scalar = core.math.Scalar;
pub const Vec2 = core.math.Vec2;
pub const Vec3 = core.math.Vec3;
pub const bond_length = core.math.bond_length;
pub const epsilon: Scalar = 0.0001;

pub fn add(a: Vec2, b: Vec2) Vec2 {
    return .{ .x = a.x + b.x, .y = a.y + b.y };
}

pub fn subtract(a: Vec2, b: Vec2) Vec2 {
    return .{ .x = a.x - b.x, .y = a.y - b.y };
}

pub fn negate(vector: Vec2) Vec2 {
    return .{ .x = -vector.x, .y = -vector.y };
}

pub fn scale(vector: Vec2, factor: Scalar) Vec2 {
    return .{ .x = vector.x * factor, .y = vector.y * factor };
}

pub fn divide(vector: Vec2, divisor: Scalar) Vec2 {
    return .{ .x = vector.x / divisor, .y = vector.y / divisor };
}

pub fn dot(a: Vec2, b: Vec2) Scalar {
    return a.x * b.x + a.y * b.y;
}

pub fn cross(a: Vec2, b: Vec2) Scalar {
    return a.x * b.y - a.y * b.x;
}

pub fn squaredLength(vector: Vec2) Scalar {
    return vector.x * vector.x + vector.y * vector.y;
}

/// Matches sketcherMinimizerPointF::length: vectors whose squared length is
/// at most SKETCHER_EPSILON have zero length.
pub fn length(vector: Vec2) Scalar {
    const squared = squaredLength(vector);
    return if (squared > epsilon) @sqrt(squared) else 0;
}

/// Matches sketcherMinimizerPointF::normalize: a degenerate vector is returned
/// unchanged rather than replaced by zero or reported as an error.
pub fn normalize(vector: Vec2) Vec2 {
    const magnitude = length(vector);
    return if (magnitude > epsilon) divide(vector, magnitude) else vector;
}

/// Rotate using upstream's screen-coordinate convention. Positive sine turns
/// (1, 0) toward (0, -1), clockwise in conventional Cartesian coordinates.
pub fn rotate(vector: Vec2, sine: Scalar, cosine: Scalar) Vec2 {
    return .{
        .x = vector.x * cosine + vector.y * sine,
        .y = -vector.x * sine + vector.y * cosine,
    };
}

/// Upstream deliberately does not guard a zero-length axis here.
pub fn parallelComponent(vector: Vec2, axis: Vec2) Vec2 {
    return divide(scale(axis, dot(vector, axis)), squaredLength(axis));
}

pub fn squaredDistance(a: Vec2, b: Vec2) Scalar {
    const dx = a.x - b.x;
    const dy = a.y - b.y;
    return dx * dx + dy * dy;
}

pub fn pointsCoincide(a: Vec2, b: Vec2) bool {
    return squaredDistance(a, b) < epsilon * epsilon;
}

/// Return whether both points are strictly within the same open half-plane.
/// A point on the line, or a degenerate line, returns false as upstream does.
pub fn sameSide(a: Vec2, b: Vec2, line_start: Vec2, line_end: Vec2) bool {
    const x = line_end.x - line_start.x;
    const y = line_end.y - line_start.y;
    if (@abs(x) > @abs(y)) {
        const slope = y / x;
        const distance_a = a.y - line_start.y - slope * (a.x - line_start.x);
        const distance_b = b.y - line_start.y - slope * (b.x - line_start.x);
        return distance_b * distance_a > 0;
    }

    const slope = x / y;
    const distance_a = a.x - line_start.x - slope * (a.y - line_start.y);
    const distance_b = b.x - line_start.x - slope * (b.y - line_start.y);
    return distance_b * distance_a > 0;
}

fn radiansToDegrees(radians: Scalar) Scalar {
    // C++ selects the float transcendental overload and rounds the
    // multiplication by 180 to float, then promotes for division by M_PI.
    const scaled: Scalar = radians * 180.0;
    return @floatCast(@as(f64, scaled) / @as(f64, std.math.pi));
}

/// Signed p1-p2-p3 angle in degrees.
pub fn signedAngle(p1: Vec2, p2: Vec2, p3: Vec2) Scalar {
    const first = subtract(p1, p2);
    const second = subtract(p3, p2);
    return radiansToDegrees(std.math.atan2(cross(first, second), dot(first, second)));
}

/// Unsigned p1-p2-p3 angle in degrees, including upstream's epsilon floor for
/// a degenerate arm and clamp against acos domain drift.
pub fn unsignedAngle(p1: Vec2, p2: Vec2, p3: Vec2) Scalar {
    const first = subtract(p1, p2);
    const second = subtract(p3, p2);
    var denominator = @sqrt(squaredLength(first)) * @sqrt(squaredLength(second));
    if (denominator < epsilon) denominator = epsilon;
    var cosine = dot(first, second) / denominator;
    if (cosine < -1) {
        cosine = -1;
    } else if (cosine > 1) {
        cosine = 1;
    }
    // Zig's software f32 acos rounds acos(-1) one ULP below the pinned C++
    // libm result; preserve upstream's exact straight-angle result.
    if (cosine == -1) return 180;
    return radiansToDegrees(std.math.acos(cosine));
}

/// Projection onto the infinite line through line_start and line_end. The
/// parameter is intentionally not clamped to the segment.
pub fn projectPointOnLine(point: Vec2, line_start: Vec2, line_end: Vec2) Vec2 {
    const from_start = subtract(point, line_start);
    const direction = subtract(line_end, line_start);
    var squared_segment_length = squaredLength(direction);
    if (squared_segment_length < epsilon) squared_segment_length = epsilon;
    const parameter = dot(from_start, direction) / squared_segment_length;
    return add(line_start, scale(direction, parameter));
}

pub const PointSegmentDistance = struct {
    squared_distance: Scalar,
    parameter: Scalar,
};

/// Squared distance to a segment with the upstream epsilon floor. `parameter`
/// is always clamped to [0, 1], unlike projectPointOnLine's parameter.
pub fn squaredDistancePointSegment(point: Vec2, start: Vec2, end: Vec2) PointSegmentDistance {
    const from_start = subtract(point, start);
    const to_end = subtract(end, point);
    const segment = subtract(end, start);
    var squared_segment_length = squaredLength(segment);
    if (squared_segment_length < epsilon) squared_segment_length = epsilon;

    const raw_parameter = dot(from_start, segment) / squared_segment_length;
    const parameter = if (raw_parameter < 0)
        @as(Scalar, 0)
    else if (raw_parameter > 1)
        @as(Scalar, 1)
    else
        raw_parameter;
    var squared_distance: Scalar = undefined;
    if (raw_parameter < 0) {
        squared_distance = squaredLength(from_start);
    } else if (raw_parameter > 1) {
        squared_distance = squaredLength(to_end);
    } else {
        squared_distance = squaredLength(subtract(point, add(start, scale(segment, raw_parameter))));
    }
    if (squared_distance < epsilon) squared_distance = epsilon;
    return .{ .squared_distance = squared_distance, .parameter = parameter };
}

/// Return the unique segment intersection. Parallel and collinear segments do
/// not intersect under the pinned upstream contract; shared endpoints do.
pub fn segmentIntersection(first_start: Vec2, first_end: Vec2, second_start: Vec2, second_end: Vec2) ?Vec2 {
    const first_direction = subtract(first_end, first_start);
    const second_direction = subtract(second_end, second_start);
    const determinant = cross(first_direction, second_direction);
    if (determinant > -epsilon and determinant < epsilon) return null;

    const between_starts = subtract(second_start, first_start);
    const first_parameter = cross(between_starts, second_direction) / determinant;
    if (first_parameter < 0 or first_parameter > 1) return null;
    const second_parameter = cross(between_starts, first_direction) / determinant;
    if (second_parameter < 0 or second_parameter > 1) return null;
    return add(first_start, scale(first_direction, first_parameter));
}

pub const Bounds2 = extern struct {
    min: Vec2,
    max: Vec2,

    pub fn isFinite(self: Bounds2) bool {
        return self.min.isFinite() and self.max.isFinite();
    }
};

/// Upstream returns a zero box for an empty molecule and otherwise seeds the
/// extrema from the first point before comparing every coordinate.
pub fn bounds(points: []const Vec2) Bounds2 {
    if (points.len == 0) return .{ .min = .{}, .max = .{} };
    var result: Bounds2 = .{ .min = points[0], .max = points[0] };
    for (points) |point| {
        if (point.x < result.min.x) result.min.x = point.x;
        if (point.y < result.min.y) result.min.y = point.y;
        if (point.x > result.max.x) result.max.x = point.x;
        if (point.y > result.max.y) result.max.y = point.y;
    }
    return result;
}

/// Arithmetic center of points, matching sketcherMinimizerMolecule::center.
pub fn center(points: []const Vec2) Vec2 {
    if (points.len == 0) return .{};
    var result: Vec2 = .{};
    for (points) |point| result = add(result, point);
    return divide(result, @floatFromInt(points.len));
}

pub const RigidTransform2 = extern struct {
    cosine: Scalar = 1,
    sine: Scalar = 0,
    translation: Vec2 = .{},

    pub fn apply(self: RigidTransform2, point: Vec2) Vec2 {
        return add(rotate(point, self.sine, self.cosine), self.translation);
    }
};

test "geometry interface conserves f32 layout" {
    try std.testing.expectEqual(@as(usize, 16), @sizeOf(Bounds2));
    try std.testing.expectEqual(@as(usize, 16), @sizeOf(RigidTransform2));
}

fn expectVec2(expected: Vec2, actual: Vec2) !void {
    try std.testing.expectApproxEqAbs(expected.x, actual.x, 0.00001);
    try std.testing.expectApproxEqAbs(expected.y, actual.y, 0.00001);
}

test "vector arithmetic preserves upstream degeneracy and rotation semantics" {
    try std.testing.expectEqual(Vec2{ .x = 4, .y = 2 }, add(.{ .x = 1, .y = 4 }, .{ .x = 3, .y = -2 }));
    try std.testing.expectEqual(Vec2{ .x = -2, .y = 6 }, subtract(.{ .x = 1, .y = 4 }, .{ .x = 3, .y = -2 }));
    try std.testing.expectEqual(@as(Scalar, 5), dot(.{ .x = 1, .y = 2 }, .{ .x = 3, .y = 1 }));
    try std.testing.expectEqual(@as(Scalar, -5), cross(.{ .x = 1, .y = 2 }, .{ .x = 3, .y = 1 }));

    const degenerate = Vec2{ .x = 0.005, .y = 0 };
    try std.testing.expectEqual(@as(Scalar, 0), length(degenerate));
    try std.testing.expectEqual(degenerate, normalize(degenerate));
    try std.testing.expectEqual(@as(Scalar, 0), length(.{ .x = 0.01 }));
    try std.testing.expect(length(.{ .x = 0.011 }) > 0);
    try expectVec2(.{ .x = 0.6, .y = 0.8 }, normalize(.{ .x = 3, .y = 4 }));
    try expectVec2(.{ .x = 0, .y = -1 }, rotate(.{ .x = 1, .y = 0 }, 1, 0));
    try expectVec2(.{ .x = 2, .y = 0 }, parallelComponent(.{ .x = 2, .y = 3 }, .{ .x = 1, .y = 0 }));
    try std.testing.expect(!parallelComponent(.{ .x = 1 }, .{}).isFinite());
}

test "angles coincidence and half-plane predicates match upstream boundaries" {
    const origin: Vec2 = .{};
    try std.testing.expectApproxEqAbs(@as(Scalar, 90), signedAngle(.{ .x = 1 }, origin, .{ .y = 1 }), 0.00001);
    try std.testing.expectApproxEqAbs(@as(Scalar, -90), signedAngle(.{ .x = 1 }, origin, .{ .y = -1 }), 0.00001);
    try std.testing.expectApproxEqAbs(@as(Scalar, 90), unsignedAngle(.{ .x = 1 }, origin, .{ .y = 1 }), 0.00001);
    try std.testing.expectEqual(@as(Scalar, 180), unsignedAngle(.{ .x = 1 }, origin, .{ .x = -1 }));
    try std.testing.expectEqual(@as(Scalar, @bitCast(@as(u32, 0x427dbd64))), signedAngle(.{ .x = 1, .y = 2 }, origin, .{ .x = -3, .y = 4 }));
    try std.testing.expectEqual(@as(Scalar, @bitCast(@as(u32, 0x427dbd64))), unsignedAngle(.{ .x = 1, .y = 2 }, origin, .{ .x = -3, .y = 4 }));
    try std.testing.expectEqual(@as(Scalar, @bitCast(@as(u32, 0xc1a84cd3))), signedAngle(.{ .x = -10, .y = -10 }, origin, .{ .x = -9, .y = -4 }));
    try std.testing.expectApproxEqAbs(@as(Scalar, @bitCast(@as(u32, 0x428f214f))), unsignedAngle(.{ .x = -10, .y = -10 }, origin, .{ .x = -8, .y = 4 }), 0.00001);
    try std.testing.expectApproxEqAbs(@as(Scalar, 90), unsignedAngle(origin, origin, .{ .x = 1 }), 0.00001);
    try std.testing.expect(std.math.isNan(unsignedAngle(.{ .x = std.math.nan(Scalar) }, origin, .{ .x = 1 })));

    try std.testing.expect(pointsCoincide(origin, .{ .x = epsilon * 0.5 }));
    try std.testing.expect(!pointsCoincide(origin, .{ .x = epsilon }));
    try std.testing.expect(sameSide(.{ .x = -1, .y = 1 }, .{ .x = 2, .y = 3 }, origin, .{ .x = 1 }));
    try std.testing.expect(!sameSide(.{ .x = -1, .y = 1 }, .{ .x = 2, .y = -3 }, origin, .{ .x = 1 }));
    try std.testing.expect(!sameSide(origin, .{ .x = 2, .y = 3 }, origin, .{ .x = 1 }));
    try std.testing.expect(sameSide(.{ .x = 1, .y = -1 }, .{ .x = 3, .y = 4 }, origin, .{ .y = 2 }));
    try std.testing.expect(!sameSide(.{ .x = 1 }, .{ .y = 1 }, origin, origin));
}

test "projection and point segment distance preserve epsilon floors and clamping" {
    try expectVec2(.{ .x = 2, .y = 0 }, projectPointOnLine(.{ .x = 2, .y = 3 }, .{}, .{ .x = 1 }));
    try expectVec2(.{ .x = 4, .y = 5 }, projectPointOnLine(.{ .x = 9, .y = 9 }, .{ .x = 4, .y = 5 }, .{ .x = 4, .y = 5 }));
    try expectVec2(.{ .x = 0.00125 }, projectPointOnLine(.{ .x = 0.005 }, .{}, .{ .x = 0.005 }));

    const middle = squaredDistancePointSegment(.{ .x = 0.5, .y = 1 }, .{}, .{ .x = 1 });
    try std.testing.expectEqual(@as(Scalar, 1), middle.squared_distance);
    try std.testing.expectEqual(@as(Scalar, 0.5), middle.parameter);
    const before = squaredDistancePointSegment(.{ .x = -2 }, .{}, .{ .x = 1 });
    try std.testing.expectEqual(@as(Scalar, 4), before.squared_distance);
    try std.testing.expectEqual(@as(Scalar, 0), before.parameter);
    const after = squaredDistancePointSegment(.{ .x = 3 }, .{}, .{ .x = 1 });
    try std.testing.expectEqual(@as(Scalar, 4), after.squared_distance);
    try std.testing.expectEqual(@as(Scalar, 1), after.parameter);
    try std.testing.expectEqual(epsilon, squaredDistancePointSegment(.{ .x = 0.5 }, .{}, .{ .x = 1 }).squared_distance);
    const nan_distance = squaredDistancePointSegment(.{ .x = std.math.nan(Scalar) }, .{}, .{ .x = 1 });
    try std.testing.expect(std.math.isNan(nan_distance.parameter));
    try std.testing.expect(std.math.isNan(nan_distance.squared_distance));
}

test "segment intersection includes endpoints and rejects parallel lines" {
    try expectVec2(.{ .x = 1, .y = 1 }, segmentIntersection(.{}, .{ .x = 2, .y = 2 }, .{ .y = 2 }, .{ .x = 2 }).?);
    try expectVec2(.{ .x = 1 }, segmentIntersection(.{}, .{ .x = 1 }, .{ .x = 1 }, .{ .x = 1, .y = 2 }).?);
    try std.testing.expect(segmentIntersection(.{}, .{ .x = 2 }, .{ .y = 1 }, .{ .x = 2, .y = 1 }) == null);
    try std.testing.expect(segmentIntersection(.{}, .{ .x = 2 }, .{ .x = 1 }, .{ .x = 3 }) == null);
    try std.testing.expect(segmentIntersection(.{}, .{ .x = 1 }, .{ .x = 2, .y = -1 }, .{ .x = 2, .y = 1 }) == null);
    try expectVec2(.{}, segmentIntersection(.{}, .{ .x = 1 }, .{}, .{ .x = 1, .y = epsilon }).?);
    try std.testing.expect(segmentIntersection(.{}, .{ .x = 1 }, .{}, .{ .x = 1, .y = epsilon * 0.5 }) == null);
}

test "bounds point center and rigid transforms retain geometric invariants" {
    const points = [_]Vec2{
        .{ .x = 2, .y = -4 },
        .{ .x = -3, .y = 8 },
        .{ .x = 7, .y = 2 },
    };
    const box = bounds(&points);
    try std.testing.expectEqual(Vec2{ .x = -3, .y = -4 }, box.min);
    try std.testing.expectEqual(Vec2{ .x = 7, .y = 8 }, box.max);
    try expectVec2(.{ .x = 2, .y = 2 }, center(&points));
    const asymmetric_points = [_]Vec2{ .{}, .{}, .{ .x = 9, .y = 3 } };
    try std.testing.expectEqual(Vec2{ .x = 3, .y = 1 }, center(&asymmetric_points));
    try std.testing.expectEqual(Bounds2{ .min = .{}, .max = .{} }, bounds(&.{}));
    try std.testing.expectEqual(Vec2{}, center(&.{}));

    const transform = RigidTransform2{ .cosine = 0, .sine = 1, .translation = .{ .x = 8, .y = -3 } };
    const first = Vec2{ .x = -2, .y = 4 };
    const second = Vec2{ .x = 5, .y = -7 };
    try std.testing.expectEqual(Vec2{ .x = 12, .y = -1 }, transform.apply(first));
    try std.testing.expectApproxEqAbs(squaredDistance(first, second), squaredDistance(transform.apply(first), transform.apply(second)), 0.00001);
    try std.testing.expectApproxEqAbs(@as(Scalar, 1), squaredLength(normalize(.{ .x = -5, .y = 12 })), 0.00001);
    try std.testing.expectEqual(squaredDistance(first, second), squaredDistance(second, first));
}
