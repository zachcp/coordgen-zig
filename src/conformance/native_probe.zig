//! The native side of the probe surface (cgz-7v2.4.2.1).
//!
//! The epic requires two surfaces around each implementation: a narrow stable
//! one, and a conformance-only one that reports internal state and is never
//! installed. `coordgen_probe.h` is the oracle's; this is native's, and it
//! lives in the conformance layer for the same reason — no production layer may
//! import this module, so nothing installed can reach it.
//!
//! It reports one seam today: the pose as it stands immediately before global
//! orientation. That is the stage the oracle's hook already captures, and
//! comparing there is what separates a layout divergence from an orientation
//! one. The remaining seams named in the parity manifest (Morgan ranks, rings,
//! fragments, DOFs, template mappings, components) belong here as they are
//! built.

const std = @import("std");
const core = @import("core");
const api = @import("api");
const generator = @import("generator");

pub const Stages = struct {
    allocator: std.mem.Allocator,
    /// Final coordinates, caller input order. Identical to what `api.generate`
    /// returns for the same input; the differential runner asserts that rather
    /// than assuming it, which makes any divergence between the two a reported
    /// nondeterminism finding instead of a silent one.
    coordinates: []core.math.Vec2,
    /// The same atoms before `orientComponents` runs.
    pre_orientation: []core.math.Vec2,
    clean_pose: bool,

    pub fn deinit(self: *Stages) void {
        self.allocator.free(self.pre_orientation);
        self.allocator.free(self.coordinates);
        self.* = undefined;
    }
};

/// Run one generation and report both stages. Validation and the pipeline are
/// the public entry point's, reached through the same `generateInto`, so this
/// cannot drift into measuring a different code path than production runs.
pub fn generateStages(allocator: std.mem.Allocator, input: api.Input) core.errors.Error!Stages {
    try input.validate();

    const coordinates = allocator.alloc(core.math.Vec2, input.atoms.len) catch
        return error.OutOfMemory;
    errdefer allocator.free(coordinates);
    const pre_orientation = allocator.alloc(core.math.Vec2, input.atoms.len) catch
        return error.OutOfMemory;
    errdefer allocator.free(pre_orientation);

    const input_to_internal = allocator.alloc(u32, input.atoms.len) catch
        return error.OutOfMemory;
    defer allocator.free(input_to_internal);
    const internal_to_input = allocator.alloc(core.ids.InputIndex, input.atoms.len) catch
        return error.OutOfMemory;
    defer allocator.free(internal_to_input);

    const clean_pose = try generator.generateInto(allocator, input, .{
        .coordinates = coordinates,
        .input_to_internal = input_to_internal,
        .internal_to_input = internal_to_input,
        .pre_orientation = pre_orientation,
    });

    return .{
        .allocator = allocator,
        .coordinates = coordinates,
        .pre_orientation = pre_orientation,
        .clean_pose = clean_pose,
    };
}

test "the probe reports a pre-orientation pose that global orientation then moves" {
    const atoms = [_]api.AtomInput{
        .{ .atomic_number = .carbon },
        .{ .atomic_number = .carbon },
        .{ .atomic_number = .carbon },
        .{ .atomic_number = .carbon },
        .{ .atomic_number = .carbon },
        .{ .atomic_number = .carbon },
        .{ .atomic_number = .oxygen },
    };
    const bonds = [_]api.BondInput{
        .{ .start = 0, .end = 1, .order = .single },
        .{ .start = 1, .end = 2, .order = .single },
        .{ .start = 2, .end = 3, .order = .single },
        .{ .start = 3, .end = 4, .order = .single },
        .{ .start = 4, .end = 5, .order = .single },
        .{ .start = 5, .end = 0, .order = .single },
        .{ .start = 2, .end = 6, .order = .single },
    };
    const input = api.Input{ .atoms = &atoms, .bonds = &bonds };

    var stages = try generateStages(std.testing.allocator, input);
    defer stages.deinit();

    try std.testing.expectEqual(atoms.len, stages.pre_orientation.len);
    for (stages.pre_orientation) |point| try std.testing.expect(point.isFinite());
    for (stages.coordinates) |point| try std.testing.expect(point.isFinite());

    // Global orientation is rigid — rotation, reflection and translation only —
    // so it cannot change any distance within a component. Comparing the two
    // stages pairwise is therefore both the strongest available invariant and
    // the assertion that separates a real captured pose from zeroed storage,
    // which would pass a finiteness check.
    for (0..atoms.len) |i| {
        for (i + 1..atoms.len) |j| {
            try std.testing.expectApproxEqAbs(
                separation(stages.pre_orientation, i, j),
                separation(stages.coordinates, i, j),
                0.01,
            );
        }
    }
    // The captured pose is a real layout, not the origin: bonded atoms sit
    // about one bond length apart in it.
    for (bonds) |bond| {
        const length = separation(stages.pre_orientation, bond.start, bond.end);
        try std.testing.expect(length > core.math.bond_length * 0.5);
        try std.testing.expect(length < core.math.bond_length * 1.5);
    }

    // And the final pose agrees with the public entry point's, which is the
    // invariant that lets the runner trust one generation for both.
    var public = try api.generate(std.testing.allocator, input);
    defer public.deinit();
    try std.testing.expectEqual(public.clean_pose, stages.clean_pose);
    for (public.coordinates, stages.coordinates) |expected, actual| {
        try std.testing.expectEqual(expected.x, actual.x);
        try std.testing.expectEqual(expected.y, actual.y);
    }
}

fn separation(coordinates: []const core.math.Vec2, left: usize, right: usize) f32 {
    const delta_x = coordinates[left].x - coordinates[right].x;
    const delta_y = coordinates[left].y - coordinates[right].y;
    return @sqrt(delta_x * delta_x + delta_y * delta_y);
}
