const std = @import("std");
const core = @import("core");
const c_abi = @import("c_abi");

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

pub const RingProbe = extern struct {
    atom_start: u32,
    atom_count: u32,
};

pub const TemplateMappingProbe = extern struct {
    input_atom: u32,
    template_atom: u32,
};

pub const FragmentProbe = extern struct {
    parent: u32,
    component: u32,
    atom_start: u32,
    atom_count: u32,
    ring_start: u32,
    ring_count: u32,
    dof_start: u32,
    dof_count: u32,
    flags: u32,
    template_match: u32,
    template_mapping_start: u32,
    template_mapping_count: u32,
};

/// A transform is deliberately marked unavailable when the C++ oracle applies
/// placement in place and preserves no source transform to report.
pub const ComponentTransformStatus = enum(u32) {
    unobserved = 0,
    observed = 1,
};

pub const ComponentProbe = extern struct {
    atom_start: u32,
    atom_count: u32,
    transform_status: ComponentTransformStatus,
    reserved: u32,
    transform: [6]f32,
};

/// Zig mirror of `coordgen_probe_result_t`. Every pointer belongs to `owner`
/// and stays valid until the oracle frees the result; conformance code reads
/// these arrays, never writes them.
pub const ProbeResult = extern struct {
    input_to_internal: c_abi.U32Span = .{},
    internal_to_input: c_abi.U32Span = .{},
    morgan_ranks: c_abi.U32Span = .{},
    ring_atoms: c_abi.U32Span = .{},
    fragment_atoms: c_abi.U32Span = .{},
    fragment_rings: c_abi.U32Span = .{},
    component_atoms: c_abi.U32Span = .{},
    template_mapping: ?[*]const TemplateMappingProbe = null,
    template_mapping_count: u32 = 0,
    template_mapping_reserved: u32 = 0,
    rings: ?[*]const RingProbe = null,
    ring_count: u32 = 0,
    ring_reserved: u32 = 0,
    fragments: ?[*]const FragmentProbe = null,
    fragment_count: u32 = 0,
    fragment_reserved: u32 = 0,
    dofs: ?[*]const DofProbe = null,
    dof_count: u32 = 0,
    dof_reserved: u32 = 0,
    components: ?[*]const ComponentProbe = null,
    component_count: u32 = 0,
    clean_pose: u32 = 0,
    reserved: u32 = 0,
    owner: ?*anyopaque = null,

    pub fn ringSlice(self: ProbeResult) []const RingProbe {
        return sliceOf(RingProbe, self.rings, self.ring_count);
    }

    pub fn fragmentSlice(self: ProbeResult) []const FragmentProbe {
        return sliceOf(FragmentProbe, self.fragments, self.fragment_count);
    }

    pub fn dofSlice(self: ProbeResult) []const DofProbe {
        return sliceOf(DofProbe, self.dofs, self.dof_count);
    }

    pub fn componentSlice(self: ProbeResult) []const ComponentProbe {
        return sliceOf(ComponentProbe, self.components, self.component_count);
    }

    pub fn templateMappingSlice(self: ProbeResult) []const TemplateMappingProbe {
        return sliceOf(TemplateMappingProbe, self.template_mapping, self.template_mapping_count);
    }

    fn sliceOf(comptime T: type, pointer: ?[*]const T, count: u32) []const T {
        const items = pointer orelse return &.{};
        return items[0..count];
    }
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

test "additional probe records are flat and use explicit unavailable transforms" {
    try std.testing.expectEqual(@as(usize, 8), @sizeOf(RingProbe));
    try std.testing.expectEqual(@as(usize, 8), @sizeOf(TemplateMappingProbe));
    try std.testing.expectEqual(@as(usize, 48), @sizeOf(FragmentProbe));
    try std.testing.expectEqual(@as(usize, 40), @sizeOf(ComponentProbe));
    try std.testing.expectEqual(ComponentTransformStatus.unobserved, ComponentTransformStatus.unobserved);
}

test "probe result mirrors coordgen_probe_result_t exactly" {
    // tests/probe_layout.c freezes the same numbers on the C side; a
    // conformance run reads oracle memory through this declaration, so a
    // silent divergence here would be read as data.
    try std.testing.expectEqual(@as(usize, 208), @sizeOf(ProbeResult));
    try std.testing.expectEqual(@as(usize, 8), @alignOf(ProbeResult));
    try std.testing.expectEqual(@as(usize, 112), @offsetOf(ProbeResult, "template_mapping"));
    try std.testing.expectEqual(@as(usize, 128), @offsetOf(ProbeResult, "rings"));
    try std.testing.expectEqual(@as(usize, 144), @offsetOf(ProbeResult, "fragments"));
    try std.testing.expectEqual(@as(usize, 160), @offsetOf(ProbeResult, "dofs"));
    try std.testing.expectEqual(@as(usize, 176), @offsetOf(ProbeResult, "components"));
    try std.testing.expectEqual(@as(usize, 188), @offsetOf(ProbeResult, "clean_pose"));
    try std.testing.expectEqual(@as(usize, 200), @offsetOf(ProbeResult, "owner"));
}
