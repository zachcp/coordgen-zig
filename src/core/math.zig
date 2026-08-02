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

test "conserved geometry representations" {
    try std.testing.expectEqual(@as(usize, 4), @sizeOf(Scalar));
    try std.testing.expectEqual(@as(usize, 8), @sizeOf(Vec2));
    try std.testing.expectEqual(@as(usize, 12), @sizeOf(Vec3));
    try std.testing.expectEqual(@as(f32, 50), bond_length);
    try std.testing.expect((Vec2{ .x = 1, .y = -2 }).isFinite());
}
