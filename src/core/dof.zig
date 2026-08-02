const std = @import("std");
const ids = @import("ids.zig");

pub const DofKind = enum(u32) {
    flip_fragment = 0,
    change_parent_bond_length = 1,
    rotate_fragment = 2,
    scale_atoms = 3,
    scale_fragment = 4,
    invert_bond = 5,
    flip_ring = 6,
};

pub const Tier = struct {
    pub const flip_fragment: u32 = 0;
    pub const invert_bond: u32 = 1;
    pub const flip_ring: u32 = 1;
    pub const change_parent_bond_length: u32 = 2;
    pub const rotate_fragment: u32 = 3;
    pub const scale_atoms: u32 = 4;
    pub const scale_fragment: u32 = 5;
};

pub const AtomRange = extern struct {
    start: u32 = 0,
    len: u32 = 0,
};

pub const State = extern struct {
    current: u32 = 0,
    optimal: u32 = 0,
    count: u32,
    tier: u32,

    pub fn validate(self: State) !void {
        if (self.count == 0 or self.current >= self.count or self.optimal >= self.count) {
            return error.InvalidDofState;
        }
    }

    pub fn advance(self: *State) void {
        std.debug.assert(self.count != 0);
        self.current = (self.current + 1) % self.count;
    }

    pub fn storeOptimal(self: *State) void {
        self.optimal = self.current;
    }

    pub fn restoreOptimal(self: *State) void {
        self.current = self.optimal;
    }
};

pub const FlipFragment = struct {};
pub const ChangeParentBondLength = struct {};
pub const RotateFragment = struct {};
pub const ScaleAtoms = struct { pivot: ids.AtomId };
pub const ScaleFragment = struct {};
pub const InvertBond = struct {
    pivot: ids.AtomId,
    bound: ids.AtomId,
};
pub const FlipRing = struct {
    ring: ids.RingId,
    pivot_a: ids.AtomId,
    pivot_b: ids.AtomId,
    /// abs(ring_size - 2 * fusion_atom_count + 2) in pinned upstream.
    penalty_multiplier: i32,
};

/// The payload is a value. Variable atom membership is represented by a range
/// into `Collection.affected_atoms`; no DOF owns pointers or allocator state.
pub const Payload = union(DofKind) {
    flip_fragment: FlipFragment,
    change_parent_bond_length: ChangeParentBondLength,
    rotate_fragment: RotateFragment,
    scale_atoms: ScaleAtoms,
    scale_fragment: ScaleFragment,
    invert_bond: InvertBond,
    flip_ring: FlipRing,
};

pub const Dof = struct {
    id: ids.DofId,
    fragment: ids.FragmentId,
    affected_atoms: AtomRange = .{},
    state: State,
    payload: Payload,

    pub fn kind(self: Dof) DofKind {
        return std.meta.activeTag(self.payload);
    }
};

/// Owned for the layout/discrete-search phases of one generation context.
/// Items are ordered by FragmentId, then pinned-upstream primary order
/// (flip/change-length/rotate), then deterministic builder discovery order.
pub const Collection = struct {
    allocator: std.mem.Allocator,
    items: []Dof,
    affected_atoms: []ids.AtomId,

    pub fn deinit(self: *Collection) void {
        self.allocator.free(self.affected_atoms);
        self.allocator.free(self.items);
        self.* = undefined;
    }
};

test "all seven pinned upstream DOF subclasses have value variants" {
    const expected = [_]DofKind{
        .flip_fragment,
        .change_parent_bond_length,
        .rotate_fragment,
        .scale_atoms,
        .scale_fragment,
        .invert_bond,
        .flip_ring,
    };
    try std.testing.expectEqual(expected.len, @typeInfo(DofKind).@"enum".field_names.len);
}

test "DOF state transitions are deterministic" {
    var state = State{ .count = 2, .tier = Tier.invert_bond };
    try state.validate();
    state.advance();
    state.storeOptimal();
    state.current = 0;
    state.restoreOptimal();
    try std.testing.expectEqual(@as(u32, 1), state.current);
}
