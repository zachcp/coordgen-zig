const std = @import("std");
const core = @import("core");

pub const Scalar = core.math.Scalar;
pub const Vec2 = core.math.Vec2;
pub const Vec3 = core.math.Vec3;
pub const bond_length = core.math.bond_length;

pub const Bounds2 = extern struct {
    min: Vec2,
    max: Vec2,

    pub fn isFinite(self: Bounds2) bool {
        return self.min.isFinite() and self.max.isFinite();
    }
};

pub const RigidTransform2 = extern struct {
    cosine: Scalar = 1,
    sine: Scalar = 0,
    translation: Vec2 = .{},
};

test "geometry interface conserves f32 layout" {
    try std.testing.expectEqual(@as(usize, 16), @sizeOf(Bounds2));
    try std.testing.expectEqual(@as(usize, 16), @sizeOf(RigidTransform2));
}
