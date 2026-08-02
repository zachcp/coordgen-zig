const std = @import("std");
const ids = @import("ids.zig");
const math = @import("math.zig");

pub const InteractionKind = enum(u32) {
    stretch = 0,
    bend = 1,
    clash = 2,
    constraint = 3,
    ez_constraint = 4,
};

pub const Stretch = struct {
    atom_a: ids.AtomId,
    atom_b: ids.AtomId,
    force_constant: f32 = 1,
    rest_length: f32 = math.bond_length,
    tolerance: f32 = 0,
};

pub const Bend = struct {
    atom_a: ids.AtomId,
    center: ids.AtomId,
    atom_b: ids.AtomId,
    force_constant: f32 = 1,
    secondary_force_constant: f32 = 0.05,
    rest_degrees: f32 = 120,
    is_ring: bool = false,
};

pub const Clash = struct {
    segment_start: ids.AtomId,
    point: ids.AtomId,
    segment_end: ids.AtomId,
    force_constant: f32 = 1,
    secondary_force_constant: f32 = 0.1,
    rest_squared_distance: f32 = 900,
};

pub const Constraint = struct {
    atom: ids.AtomId,
    origin: math.Vec2,
    force_constant: f32 = 0.5,
};

pub const EzConstraint = struct {
    side_a: ids.AtomId,
    double_a: ids.AtomId,
    double_b: ids.AtomId,
    side_b: ids.AtomId,
    is_z: bool,
    force_movement: bool = false,
};

/// Every concrete `sketcherMinimizerInteraction` subclass constructed by the
/// pinned source is represented as a value-owned tagged union.
pub const Payload = union(InteractionKind) {
    stretch: Stretch,
    bend: Bend,
    clash: Clash,
    constraint: Constraint,
    ez_constraint: EzConstraint,
};

pub const Interaction = struct {
    id: ids.InteractionId,
    payload: Payload,

    pub fn kind(self: Interaction) InteractionKind {
        return std.meta.activeTag(self.payload);
    }
};

pub const Collection = struct {
    allocator: std.mem.Allocator,
    items: []Interaction,

    pub fn deinit(self: *Collection) void {
        self.allocator.free(self.items);
        self.* = undefined;
    }
};

test "all pinned upstream minimizer interaction subclasses are covered" {
    try std.testing.expectEqual(@as(usize, 5), @typeInfo(InteractionKind).@"enum".field_names.len);
    const item = Interaction{
        .id = ids.InteractionId.fromIndex(0),
        .payload = .{ .stretch = .{
            .atom_a = ids.AtomId.fromIndex(0),
            .atom_b = ids.AtomId.fromIndex(1),
        } },
    };
    try std.testing.expectEqual(InteractionKind.stretch, item.kind());
}
