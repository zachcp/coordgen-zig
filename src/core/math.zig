const std = @import("std");

/// Initial parity intentionally uses f32, matching upstream evaluation.
pub const Scalar = f32;
pub const bond_length: Scalar = 50.0;

pub const Vec2 = extern struct {
    x: Scalar = 0,
    y: Scalar = 0,

    pub fn isFinite(self: Vec2) bool {
        return std.math.isFinite(self.x) and std.math.isFinite(self.y);
    }
};

pub const Vec3 = extern struct {
    x: Scalar = 0,
    y: Scalar = 0,
    z: Scalar = 0,

    pub fn isFinite(self: Vec3) bool {
        return std.math.isFinite(self.x) and
            std.math.isFinite(self.y) and
            std.math.isFinite(self.z);
    }
};

/// Upstream's `roundToTwoDecimalDigits(float)`: multiplication by 100 occurs
/// in float, while the unsuffixed literals promote the remaining expression
/// through double before the final cast back to float.
pub fn roundToTwoDecimalDigits(value: f32) f32 {
    const scaled: f32 = value * 100;
    return @floatCast(@floor(@as(f64, scaled) + 0.5) * 0.01);
}

test "conserved geometry representations" {
    try std.testing.expectEqual(@as(usize, 4), @sizeOf(Scalar));
    try std.testing.expectEqual(@as(usize, 8), @sizeOf(Vec2));
    try std.testing.expectEqual(@as(usize, 12), @sizeOf(Vec3));
    try std.testing.expectEqual(@as(f32, 50), bond_length);
    try std.testing.expect((Vec2{ .x = 1, .y = -2 }).isFinite());
}

test "hundredth rounding preserves upstream mixed precision" {
    try std.testing.expectEqual(@as(u32, 0xc2c7f5c3), @as(u32, @bitCast(roundToTwoDecimalDigits(-99.985))));
}
