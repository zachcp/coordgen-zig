const std = @import("std");
const core = @import("core");

/// Stable only within a matching oracle/native conformance build. This module
/// and its C header are never part of installed production artifacts.
pub const DofProbe = extern struct {
    id: u32,
    kind: core.dof.DofKind,
    fragment: u32,
    current_state: u32,
    optimal_state: u32,
    state_count: u32,
    tier: u32,
    affected_start: u32,
    affected_count: u32,
    atom_a: u32,
    atom_b: u32,
    ring: u32,
    current_penalty: f32,
    variant_penalty_multiplier: i32,
};

pub fn dofProbe(dof: core.dof.Dof, current_penalty: f32) DofProbe {
    var probe = DofProbe{
        .id = dof.id.index(),
        .kind = dof.kind(),
        .fragment = dof.fragment.index(),
        .current_state = dof.state.current,
        .optimal_state = dof.state.optimal,
        .state_count = dof.state.count,
        .tier = dof.state.tier,
        .affected_start = dof.affected_atoms.start,
        .affected_count = dof.affected_atoms.len,
        .atom_a = @backingInt(core.ids.AtomId.invalid),
        .atom_b = @backingInt(core.ids.AtomId.invalid),
        .ring = @backingInt(core.ids.RingId.invalid),
        .current_penalty = current_penalty,
        .variant_penalty_multiplier = 0,
    };
    switch (dof.payload) {
        .scale_atoms => |value| probe.atom_a = value.pivot.index(),
        .invert_bond => |value| {
            probe.atom_a = value.pivot.index();
            probe.atom_b = value.bound.index();
        },
        .flip_ring => |value| {
            probe.atom_a = value.pivot_a.index();
            probe.atom_b = value.pivot_b.index();
            probe.ring = value.ring.index();
            probe.variant_penalty_multiplier = value.penalty_multiplier;
        },
        else => {},
    }
    return probe;
}

test "DOF probe is flat and preserves state, penalty, IDs, and membership" {
    const value = core.dof.Dof{
        .id = core.ids.DofId.fromIndex(4),
        .fragment = core.ids.FragmentId.fromIndex(2),
        .affected_atoms = .{ .start = 8, .len = 3 },
        .state = .{ .current = 1, .optimal = 1, .count = 2, .tier = core.dof.Tier.invert_bond },
        .payload = .{ .invert_bond = .{
            .pivot = core.ids.AtomId.fromIndex(5),
            .bound = core.ids.AtomId.fromIndex(6),
        } },
    };
    const probe = dofProbe(value, 100);
    try std.testing.expectEqual(@as(usize, 56), @sizeOf(DofProbe));
    try std.testing.expectEqual(core.dof.DofKind.invert_bond, probe.kind);
    try std.testing.expectEqual(@as(u32, 5), probe.atom_a);
    try std.testing.expectEqual(@as(u32, 6), probe.atom_b);
    try std.testing.expectEqual(@as(f32, 100), probe.current_penalty);
}
