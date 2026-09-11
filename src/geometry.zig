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
    // The pinned Clang build contracts upstream's multiply-add expression.
    return @mulAdd(Scalar, a.x, b.x, a.y * b.y);
}

pub fn cross(a: Vec2, b: Vec2) Scalar {
    // Keep parity with Clang's contracted multiply-subtract in PointF.
    return @mulAdd(Scalar, a.x, b.y, -(a.y * b.x));
}

pub fn squaredLength(vector: Vec2) Scalar {
    // Keep parity with Clang's contracted multiply-add in PointF.
    return @mulAdd(Scalar, vector.x, vector.x, vector.y * vector.y);
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
    // PointF::rotate contracts each multiply-add in the pinned Clang build.
    return .{
        .x = @mulAdd(Scalar, vector.x, cosine, vector.y * sine),
        .y = @mulAdd(Scalar, -vector.x, sine, vector.y * cosine),
    };
}

/// Upstream deliberately does not guard a zero-length axis here.
pub fn parallelComponent(vector: Vec2, axis: Vec2) Vec2 {
    return divide(scale(axis, dot(vector, axis)), squaredLength(axis));
}

pub fn squaredDistance(a: Vec2, b: Vec2) Scalar {
    const dx = a.x - b.x;
    const dy = a.y - b.y;
    // PointF::squaredDistance contracts the sum in the pinned Clang build.
    return @mulAdd(Scalar, dx, dx, dy * dy);
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

// Derived from Sun Microsystems' fdlibm e_acosf.c, used by glibc 2.36.
// Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
// Developed at SunPro, a Sun Microsystems, Inc. business.
// Permission to use, copy, modify, and distribute this software is freely
// granted, provided that this notice is preserved.
fn upstreamAcos(value: Scalar) Scalar {
    const one: Scalar = 1;
    const pi: Scalar = @bitCast(@as(u32, 0x40490fda));
    const pio2_hi: Scalar = @bitCast(@as(u32, 0x3fc90fda));
    const pio2_lo: Scalar = @bitCast(@as(u32, 0x33a22168));
    const p_s0: Scalar = @bitCast(@as(u32, 0x3e2aaaab));
    const p_s1: Scalar = @bitCast(@as(u32, 0xbea6b090));
    const p_s2: Scalar = @bitCast(@as(u32, 0x3e4e0aa8));
    const p_s3: Scalar = @bitCast(@as(u32, 0xbd241146));
    const p_s4: Scalar = @bitCast(@as(u32, 0x3a4f7f04));
    const p_s5: Scalar = @bitCast(@as(u32, 0x3811ef08));
    const q_s1: Scalar = @bitCast(@as(u32, 0xc019d139));
    const q_s2: Scalar = @bitCast(@as(u32, 0x4001572d));
    const q_s3: Scalar = @bitCast(@as(u32, 0xbf303361));
    const q_s4: Scalar = @bitCast(@as(u32, 0x3d9dc62e));

    const bits: u32 = @bitCast(value);
    const magnitude = bits & 0x7fffffff;
    if (magnitude == 0x3f800000) {
        if (bits >> 31 == 0) return 0;
        return pi + 2 * pio2_lo;
    }
    if (magnitude > 0x3f800000) return std.math.nan(Scalar);

    if (magnitude < 0x3f000000) {
        if (magnitude <= 0x32800000) return pio2_hi + pio2_lo;
        const z = value * value;
        const p = z * (p_s0 + z * (p_s1 + z * (p_s2 + z * (p_s3 + z * (p_s4 + z * p_s5)))));
        const q = one + z * (q_s1 + z * (q_s2 + z * (q_s3 + z * q_s4)));
        const r = p / q;
        return pio2_hi - (value - (pio2_lo - value * r));
    }

    const z = if (bits >> 31 != 0) (one + value) * 0.5 else (one - value) * 0.5;
    const p = z * (p_s0 + z * (p_s1 + z * (p_s2 + z * (p_s3 + z * (p_s4 + z * p_s5)))));
    const q = one + z * (q_s1 + z * (q_s2 + z * (q_s3 + z * q_s4)));
    const s = @sqrt(z);
    const r = p / q;
    if (bits >> 31 != 0) {
        const w = r * s - pio2_lo;
        return pi - 2 * (s + w);
    }

    const df: Scalar = @bitCast(@as(u32, @bitCast(s)) & 0xfffff000);
    const c = (z - df * df) / (s + df);
    const w = r * s + c;
    return 2 * (df + w);
}

/// Signed p1-p2-p3 angle in degrees.
pub fn signedAngle(p1: Vec2, p2: Vec2, p3: Vec2) Scalar {
    const first = subtract(p1, p2);
    const second = subtract(p3, p2);
    // Clang contracts upstream's `x1 * y2 - y1 * x2` expression into an FMA.
    // Keep that single-rounding residue: even the angle from a vector to
    // itself can be microscopically negative, which anchors clockwise ties.
    const cross_product = @mulAdd(Scalar, first.x, second.y, -(first.y * second.x));
    return radiansToDegrees(std.math.atan2(cross_product, dot(first, second)));
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
    return radiansToDegrees(upstreamAcos(cosine));
}

/// Projection onto the infinite line through line_start and line_end. The
/// parameter is intentionally not clamped to the segment.
pub fn projectPointOnLine(point: Vec2, line_start: Vec2, line_end: Vec2) Vec2 {
    const from_start = subtract(point, line_start);
    const direction = subtract(line_end, line_start);
    var squared_segment_length = squaredLength(direction);
    if (squared_segment_length < epsilon) squared_segment_length = epsilon;
    const parameter = dot(from_start, direction) / squared_segment_length;
    // PointF's scale and add operators materialize the product before the sum;
    // unlike direct component arithmetic, this expression is not contracted.
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
        // Preserve the same operator boundary as projectPointOnLine above.
        const projection = add(start, scale(segment, raw_parameter));
        const from_projection = subtract(point, projection);
        squared_distance = squaredLength(from_projection);
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
    try std.testing.expectEqual(
        @as(Scalar, @bitCast(@as(u32, 0xc9d1a1e6))),
        dot(.{ .x = -1757.9722, .y = -700.3646 }, .{ .x = 1009.6434, .y = -82.26597 }),
    );
    try std.testing.expectEqual(@as(Scalar, -5), cross(.{ .x = 1, .y = 2 }, .{ .x = 3, .y = 1 }));
    try std.testing.expectEqual(
        @as(Scalar, @bitCast(@as(u32, 0x467dcab5))),
        cross(
            .{ .x = -154.942642, .y = 181.409866 },
            .{ .x = -181.409866, .y = 107.567947 },
        ),
    );
    try std.testing.expectEqual(
        @as(Scalar, @bitCast(@as(u32, 0x4671db1b))),
        squaredLength(.{ .x = -123.949234, .y = -10.7407618 }),
    );
    try std.testing.expectEqual(
        @as(Scalar, @bitCast(@as(u32, 0x467dcab5))),
        rotate(.{ .x = -154.942642, .y = 181.409866 }, 181.409866, 107.567947).x,
    );
    try std.testing.expectEqual(
        @as(Scalar, @bitCast(@as(u32, 0x4671db1b))),
        squaredDistance(.{ .x = -123.949234, .y = -10.7407618 }, .{}),
    );

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
    try std.testing.expectEqual(
        @as(Scalar, @bitCast(@as(u32, 0xb4954269))),
        signedAngle(.{ .x = 64.34, .y = -27.58 }, origin, .{ .x = 64.34, .y = -27.58 }),
    );
    try std.testing.expectEqual(
        @as(Scalar, @bitCast(@as(u32, 0x426fff42))),
        unsignedAngle(.{ .x = 70.05, .y = -87.84 }, .{ .x = 20.05, .y = -87.84 }, .{ .x = 45.05, .y = -131.14 }),
    );
    try std.testing.expectEqual(@as(Scalar, 0x1.0c145p0), upstreamAcos(0x1.00017p-1));
    try std.testing.expectEqual(@as(Scalar, 0x1.02a348p0), upstreamAcos(0x1.102e5p-1));
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

    const point = Vec2{ .x = -18.88, .y = -6.61 };
    const segment_start = Vec2{ .x = -75.46, .y = -21.17 };
    const segment_end = Vec2{ .x = -33.55, .y = 34.9 };
    try std.testing.expectEqual(
        Vec2{ .x = @bitCast(@as(u32, 0xc240ca3a)), .y = @bitCast(@as(u32, 0x4174db96)) },
        projectPointOnLine(point, segment_start, segment_end),
    );
    try std.testing.expectEqual(
        @as(Scalar, @bitCast(@as(u32, 0x44a77716))),
        squaredDistancePointSegment(point, segment_start, segment_end).squared_distance,
    );
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

test "vector and rigid-transform properties hold across a deterministic grid" {
    const values = [_]Scalar{ -100, -7, -0.25, 0, 0.25, 7, 100 };
    const quarter_turns = [_]RigidTransform2{
        .{ .cosine = 1, .sine = 0, .translation = .{ .x = 13, .y = -9 } },
        .{ .cosine = 0, .sine = 1, .translation = .{ .x = -4, .y = 11 } },
        .{ .cosine = -1, .sine = 0, .translation = .{ .x = 5, .y = 3 } },
        .{ .cosine = 0, .sine = -1, .translation = .{ .x = -8, .y = -2 } },
    };

    for (values) |x1| for (values) |y1| {
        const first: Vec2 = .{ .x = x1, .y = y1 };
        try std.testing.expectEqual(first, negate(negate(first)));
        if (length(first) > epsilon) {
            try std.testing.expectApproxEqAbs(@as(Scalar, 1), squaredLength(normalize(first)), 0.00001);
        }

        for (values) |x2| for (values) |y2| {
            const second: Vec2 = .{ .x = x2, .y = y2 };
            try std.testing.expectEqual(first, subtract(add(first, second), second));
            try std.testing.expectEqual(
                squaredDistance(first, second),
                squaredDistance(second, first),
            );
            try std.testing.expect(squaredDistance(first, second) >= 0);
            for (quarter_turns) |transform| {
                try std.testing.expectApproxEqAbs(
                    squaredDistance(first, second),
                    squaredDistance(transform.apply(first), transform.apply(second)),
                    0.01,
                );
            }
        };
    };
}
