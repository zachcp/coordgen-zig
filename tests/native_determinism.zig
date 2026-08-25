//! Generation-level determinism, concurrency and allocation-failure coverage
//! (cgz-7v2.4.3).
//!
//! PR #14 gave the individual seams their property and malformed-input
//! families. What it could not give them was coverage of `api.generate` as a
//! whole, because the pipeline was not wired to the public entry point yet
//! (cgz-7v2.21). This file is that half: it takes one small representative of
//! every generation path native owns and asserts four properties of the whole
//! call.
//!
//! The four are deliberately the ones that hold UNIVERSALLY, so a failure is
//! always a defect rather than a judgement call:
//!
//!   repeat            the same input twice produces the same bytes
//!   allocator-order   a different allocator does not change the result
//!   thread            concurrent independent generations agree with a serial one
//!   allocation failure  every reachable allocation index is reported and cleaned
//!
//! Input-order invariance is deliberately NOT in that list. Upstream's own
//! layout is pointer-order dependent for some inputs - that is what the parity
//! ceiling exists for (cgz-r13) - so asserting native is order-invariant in
//! general would be asserting something untrue of the thing being ported.
//!
//! What the fifth test asserts instead is the universal half - a permuted
//! input is itself deterministic - and then pins, per family, which of three
//! things reversal actually does. Measuring that turned out to be worth more
//! than a boolean would have been: of the six families that are not invariant,
//! four are merely CONGRUENT, the same shape reflected or turned, and only
//! two produce a materially different layout. A boolean would have filed all
//! six under the same heading and hidden which two matter.

const std = @import("std");
const api = @import("api");
const core = @import("core");

/// One small member per generation path. Small on purpose: the allocation
/// sweep re-runs a whole generation once per reachable allocation index, so
/// member size shows up quadratically in this file's runtime.
const Family = struct {
    name: []const u8,
    atoms: []const api.AtomInput,
    bonds: []const api.BondInput,
    residues: []const api.ResidueInput = &.{},
    residue_interactions: []const api.ResidueInteractionInput = &.{},
    /// What reversing the caller's atom order does to the layout. Measured,
    /// not assumed - see the header and `OrderStability`.
    order: OrderStability,

    fn input(self: Family) api.Input {
        return .{
            .atoms = self.atoms,
            .bonds = self.bonds,
            .residues = self.residues,
            .residue_interactions = self.residue_interactions,
        };
    }
};

/// How a layout responds to the caller reversing its atom order. Three levels
/// rather than a boolean, because "not invariant" turned out to cover two very
/// different things and collapsing them would hide the interesting one.
const OrderStability = enum {
    /// The reversed run reproduces the forward layout exactly under the
    /// reversal.
    identical,
    /// Every pairwise distance is preserved within the tolerant-tier bound, so
    /// it is the same shape placed or turned differently - a reflection or a
    /// rotation, not a different answer.
    congruent,
    /// A materially different layout. Upstream is pointer-order dependent for
    /// some inputs too, which is what the parity ceiling exists for (cgz-r13),
    /// so this is recorded rather than treated as a defect.
    divergent,
};

/// The tolerant tier's bound, 0.1 bond lengths, from
/// docs/architecture/SUCCESS_CRITERIA.md. Used here as the project's own
/// definition of "the same layout" rather than a threshold invented for this
/// file.
const congruence_tolerance: f32 = api.bond_length * 0.1;

const carbon: api.AtomInput = .{};

fn chainBonds(comptime count: u32) [count - 1]api.BondInput {
    var bonds: [count - 1]api.BondInput = undefined;
    for (&bonds, 0..) |*bond, index| bond.* = .{ .start = @intCast(index), .end = @intCast(index + 1) };
    return bonds;
}

fn ringBonds(comptime count: u32) [count]api.BondInput {
    var bonds: [count]api.BondInput = undefined;
    for (&bonds, 0..) |*bond, index| bond.* = .{
        .start = @intCast(index),
        .end = @intCast((index + 1) % count),
    };
    return bonds;
}

const chain_atoms = [_]api.AtomInput{ carbon, carbon, carbon };
const chain_bonds = chainBonds(3);

const ring_atoms = [_]api.AtomInput{ carbon, carbon, carbon, carbon, carbon, carbon };
const ring_bonds = ringBonds(6);

const fused_atoms = [_]api.AtomInput{ carbon, carbon, carbon, carbon, carbon, carbon };
const fused_bonds = [_]api.BondInput{
    .{ .start = 0, .end = 1 },
    .{ .start = 1, .end = 2 },
    .{ .start = 2, .end = 3 },
    .{ .start = 3, .end = 0 },
    .{ .start = 2, .end = 5 },
    .{ .start = 5, .end = 4 },
    .{ .start = 4, .end = 3 },
};

const macrocycle_atoms = [_]api.AtomInput{
    carbon, carbon, carbon, carbon, carbon, carbon, carbon,
    carbon, carbon, carbon, carbon, carbon, carbon,
};
const macrocycle_bonds = ringBonds(13);

/// Two chains with no bond and no relation between them: the plain component
/// arrangement path.
const component_atoms = [_]api.AtomInput{ carbon, carbon, carbon, carbon };
const component_bonds = [_]api.BondInput{
    .{ .start = 0, .end = 1 },
    .{ .start = 2, .end = 3 },
};

/// Two molecules joined by two zero-order bonds to the same partner atom. Two
/// relations rather than one is what lets the proximity score term be non-zero
/// (cgz-7v2.24), so this member exercises that branch as well as the proximity
/// arrangement.
const proximity_atoms = [_]api.AtomInput{ carbon, carbon, carbon, carbon, carbon, carbon, carbon, carbon };
const proximity_bonds = [_]api.BondInput{
    .{ .start = 0, .end = 1 },
    .{ .start = 1, .end = 2 },
    .{ .start = 2, .end = 3 },
    .{ .start = 3, .end = 4 },
    .{ .start = 4, .end = 5 },
    .{ .start = 6, .end = 7 },
    .{ .start = 0, .end = 6, .order = .zero },
    .{ .start = 5, .end = 6, .order = .zero },
};

const counterion_atoms = [_]api.AtomInput{
    .{ .atomic_number = .nitrogen, .formal_charge = 1 },
    carbon,
    .{ .atomic_number = .chlorine, .formal_charge = -1 },
};
const counterion_bonds = [_]api.BondInput{.{ .start = 0, .end = 1 }};

/// A nonterminal bond to a metal, which preparation rewrites to zero order
/// (cgz-r11). Placed on a ring so the rewrite changes the graph's shape rather
/// than only an observable.
const metal_atoms = [_]api.AtomInput{
    carbon,
    carbon,
    .{ .atomic_number = .iron },
    carbon,
    .{ .atomic_number = .oxygen },
};
const metal_bonds = [_]api.BondInput{
    .{ .start = 0, .end = 1 },
    .{ .start = 1, .end = 2 },
    .{ .start = 2, .end = 3 },
    .{ .start = 3, .end = 4 },
};

const hetero_atoms = [_]api.AtomInput{
    .{ .atomic_number = .carbon },
    .{ .atomic_number = .oxygen },
    .{ .atomic_number = .nitrogen },
    .{ .atomic_number = .sulfur },
    .{ .atomic_number = .fluorine },
};
const hetero_bonds = [_]api.BondInput{
    .{ .start = 0, .end = 1, .order = .double },
    .{ .start = 1, .end = 2 },
    .{ .start = 2, .end = 3, .order = .triple },
    .{ .start = 3, .end = 4 },
};

const residue_atoms = [_]api.AtomInput{ carbon, carbon, carbon, carbon };
const residue_bonds = [_]api.BondInput{ .{ .start = 0, .end = 1 }, .{ .start = 1, .end = 2 } };
const residue_residues = [_]api.ResidueInput{
    .{ .atom = 3, .chain = "A", .residue_number = 9, .closest_ligand_atom = 1 },
};
const residue_interactions = [_]api.ResidueInteractionInput{.{ .start = 3, .end = 1 }};

/// `order` is measured, not chosen. The worst pairwise-distance deviation under
/// reversal, in units (one bond is 50), as measured on aarch64-macos Debug when
/// this table was written:
///
///   chain             0.000000  identical
///   ring              0.000015  congruent
///   fused rings       0.000008  congruent
///   macrocycle       84.492200  divergent
///   components        0.071243  congruent
///   proximity        83.474560  divergent
///   counterion        0.000000  identical
///   metal zero-order  1.460194  congruent
///   heteroatoms       0.000000  identical
///
/// Recorded so a reclassification can be compared against a number rather than
/// argued about. `metal zero-order` is the one with the least headroom, at
/// 1.46 against a 5.0 bound.
const families = [_]Family{
    .{ .name = "chain", .atoms = &chain_atoms, .bonds = &chain_bonds, .order = .identical },
    .{ .name = "ring", .atoms = &ring_atoms, .bonds = &ring_bonds, .order = .congruent },
    .{ .name = "fused rings", .atoms = &fused_atoms, .bonds = &fused_bonds, .order = .congruent },
    .{ .name = "macrocycle", .atoms = &macrocycle_atoms, .bonds = &macrocycle_bonds, .order = .divergent },
    .{ .name = "components", .atoms = &component_atoms, .bonds = &component_bonds, .order = .congruent },
    .{ .name = "proximity", .atoms = &proximity_atoms, .bonds = &proximity_bonds, .order = .divergent },
    .{ .name = "counterion", .atoms = &counterion_atoms, .bonds = &counterion_bonds, .order = .identical },
    .{ .name = "metal zero-order", .atoms = &metal_atoms, .bonds = &metal_bonds, .order = .congruent },
    .{ .name = "heteroatoms", .atoms = &hetero_atoms, .bonds = &hetero_bonds, .order = .identical },
    .{
        .name = "residue",
        .atoms = &residue_atoms,
        .bonds = &residue_bonds,
        .residues = &residue_residues,
        .residue_interactions = &residue_interactions,
        .order = .identical,
    },
};

fn distance(first: api.Vec2, second: api.Vec2) f32 {
    const dx = first.x - second.x;
    const dy = first.y - second.y;
    return @sqrt(dx * dx + dy * dy);
}

/// Every observable the public result publishes, so "the same result" means
/// all of it and not just the coordinates.
fn expectSameResult(expected: api.Result, actual: api.Result) !void {
    try std.testing.expectEqual(expected.clean_pose, actual.clean_pose);
    try std.testing.expectEqualSlices(api.Vec2, expected.coordinates, actual.coordinates);
    try std.testing.expectEqualSlices(u32, expected.input_to_internal, actual.input_to_internal);
    try std.testing.expectEqualSlices(api.InputIndex, expected.internal_to_input, actual.internal_to_input);
    try std.testing.expectEqualSlices(api.BondOrder, expected.effective_bond_orders, actual.effective_bond_orders);
    try std.testing.expectEqualSlices(api.BondDisplay, expected.bond_displays, actual.bond_displays);
    try std.testing.expectEqualSlices(api.AtomStereo, expected.atom_stereo, actual.atom_stereo);
}

test "every generation family repeats byte for byte" {
    for (families) |family| {
        errdefer std.debug.print("family: {s}\n", .{family.name});
        var first = try api.generate(std.testing.allocator, family.input());
        defer first.deinit();
        var second = try api.generate(std.testing.allocator, family.input());
        defer second.deinit();
        try expectSameResult(first, second);
        for (first.coordinates) |coordinate| try std.testing.expect(coordinate.isFinite());
    }
}

test "every generation family is independent of the allocator it runs on" {
    for (families) |family| {
        errdefer std.debug.print("family: {s}\n", .{family.name});
        var baseline = try api.generate(std.testing.allocator, family.input());
        defer baseline.deinit();

        // A fresh arena has a different residual heap layout and hands back
        // different addresses in a different order. Nothing about the result
        // may depend on either.
        var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
        defer arena.deinit();
        const arena_result = try api.generate(arena.allocator(), family.input());
        try expectSameResult(baseline, arena_result);

        // A failing allocator that never actually fails still perturbs the
        // path taken through std's allocation bookkeeping.
        var counting = std.testing.FailingAllocator.init(std.testing.allocator, .{});
        var counted = try api.generate(counting.allocator(), family.input());
        try expectSameResult(baseline, counted);
        counted.deinit();
        // With the result released, the perturbing allocator must be square:
        // the run leaked nothing on the success path either.
        try std.testing.expectEqual(counting.allocated_bytes, counting.freed_bytes);
    }
}

const ThreadContext = struct {
    family: Family,
    /// Compared against the serial result rather than only against the other
    /// threads: agreeing with each other while both differing from a serial
    /// run would be a shared-state defect that thread-to-thread comparison
    /// cannot see.
    expected: *const api.Result,
    failure: ?anyerror = null,
    mismatch: bool = false,
};

fn generateOnThread(context: *ThreadContext) void {
    // page_allocator rather than the testing allocator: std.testing.allocator
    // is not itself the subject here, and sharing it across threads would test
    // the harness instead of the library.
    var result = api.generate(std.heap.page_allocator, context.family.input()) catch |err| {
        context.failure = err;
        return;
    };
    defer result.deinit();
    expectSameResult(context.expected.*, result) catch {
        context.mismatch = true;
    };
}

test "every generation family agrees with itself across four concurrent contexts" {
    for (families) |family| {
        errdefer std.debug.print("family: {s}\n", .{family.name});
        var serial = try api.generate(std.testing.allocator, family.input());
        defer serial.deinit();

        var contexts: [4]ThreadContext = undefined;
        var threads: [contexts.len]std.Thread = undefined;
        for (&contexts) |*context| context.* = .{ .family = family, .expected = &serial };
        for (&threads, &contexts) |*thread, *context| {
            thread.* = try std.Thread.spawn(.{}, generateOnThread, .{context});
        }
        for (&threads) |thread| thread.join();
        for (&contexts) |context| {
            if (context.failure) |err| return err;
            try std.testing.expect(!context.mismatch);
        }
    }
}

fn generateAndDiscard(allocator: std.mem.Allocator, family: Family) !void {
    var result = try api.generate(allocator, family.input());
    result.deinit();
}

test "every generation family reports and cleans up failure at every allocation index" {
    for (families) |family| {
        errdefer std.debug.print("family: {s}\n", .{family.name});
        // The counting run inside checkAllocationFailures establishes how many
        // indices exist, so a first-call warm-up that later calls skip would
        // leave the last index unexercised. Discharge it here.
        try generateAndDiscard(std.testing.allocator, family);
        // The success path must be square before failing anything means
        // anything.
        var leak_check = std.testing.FailingAllocator.init(std.testing.allocator, .{});
        try generateAndDiscard(leak_check.allocator(), family);
        try std.testing.expectEqual(leak_check.allocated_bytes, leak_check.freed_bytes);
        // No second count comparison here on purpose. A raw alloc_index taken
        // against the shared testing allocator is not a reproducible quantity -
        // whether an ArrayList grows in place depends on the residual heap, and
        // this suite measured 227 then 225 on the metal family for that reason.
        // core.oom.checkAllocationFailures counts under an allocator that
        // declines in-place growth, so its count IS stable and a count that
        // still varies is reported as NondeterministicMemoryUsage. The
        // stability assertion is kept; the unreproducible one is not
        // duplicated.
        try core.oom.checkAllocationFailures(std.testing.allocator, generateAndDiscard, .{family});
    }
}

/// The reversal remap: caller index `i` of the forward input is index
/// `len - 1 - i` of the reversed one.
fn reverseInput(
    comptime max_atoms: usize,
    comptime max_bonds: usize,
    family: Family,
    atoms: *[max_atoms]api.AtomInput,
    bonds: *[max_bonds]api.BondInput,
) api.Input {
    const count: u32 = @intCast(family.atoms.len);
    for (family.atoms, 0..) |atom, index| atoms[count - 1 - index] = atom;
    for (family.bonds, 0..) |bond, index| bonds[index] = .{
        .start = count - 1 - bond.start,
        .end = count - 1 - bond.end,
        .order = bond.order,
        .skip = bond.skip,
        .stereo = bond.stereo,
        .display = bond.display,
        .crossing_penalty_multiplier = bond.crossing_penalty_multiplier,
    };
    return .{ .atoms = atoms[0..family.atoms.len], .bonds = bonds[0..family.bonds.len] };
}

test "a reversed caller order is itself deterministic, and invariance is pinned per family" {
    const max_atoms = 16;
    const max_bonds = 16;
    for (families) |family| {
        // The reversal remap here does not carry residue references, so the
        // residue family keeps its own dedicated order test in
        // tests/native_minimal.zig.
        if (family.residues.len != 0) continue;
        errdefer std.debug.print("family: {s}\n", .{family.name});
        var atoms: [max_atoms]api.AtomInput = undefined;
        var bonds: [max_bonds]api.BondInput = undefined;
        const reversed_input = reverseInput(max_atoms, max_bonds, family, &atoms, &bonds);

        var reversed = try api.generate(std.testing.allocator, reversed_input);
        defer reversed.deinit();
        var reversed_again = try api.generate(std.testing.allocator, reversed_input);
        defer reversed_again.deinit();
        // This is the universal half: whatever a permutation produces, it
        // produces it every time.
        try expectSameResult(reversed, reversed_again);

        var forward = try api.generate(std.testing.allocator, family.input());
        defer forward.deinit();
        try std.testing.expectEqual(forward.coordinates.len, reversed.coordinates.len);

        var identical = forward.clean_pose == reversed.clean_pose;
        for (forward.coordinates, 0..) |coordinate, index| {
            const mirrored = reversed.coordinates[forward.coordinates.len - 1 - index];
            identical = identical and coordinate.x == mirrored.x and coordinate.y == mirrored.y;
        }

        // Congruence is judged on the pairwise distance matrix, which is
        // invariant under every rotation, translation and reflection at once.
        // That is what separates "turned around" from "laid out differently"
        // without having to recover the transform.
        var worst_distance_delta: f32 = 0;
        for (forward.coordinates, 0..) |first, i| {
            for (forward.coordinates, 0..) |second, j| {
                const forward_distance = distance(first, second);
                const reversed_distance = distance(
                    reversed.coordinates[forward.coordinates.len - 1 - i],
                    reversed.coordinates[forward.coordinates.len - 1 - j],
                );
                worst_distance_delta = @max(worst_distance_delta, @abs(forward_distance - reversed_distance));
            }
        }
        const congruent = forward.clean_pose == reversed.clean_pose and
            worst_distance_delta <= congruence_tolerance;

        const measured: OrderStability = if (identical)
            .identical
        else if (congruent)
            .congruent
        else
            .divergent;
        errdefer std.debug.print(
            "worst pairwise distance deviation: {d:.6} units, tolerance {d:.6}\n",
            .{ worst_distance_delta, congruence_tolerance },
        );
        // Pinned in both directions. A family that drops a level is a
        // regression; one that gains a level is a real improvement, and should
        // be recorded in the table above rather than silently absorbed.
        try std.testing.expectEqual(family.order, measured);
    }
}

/// The committed floors. Embedded rather than read from disk so the test has
/// no working-directory assumption and the table travels with the binary. The
/// parser and the comparison live in core.oom, shared with the C ABI suite.
const allocation_site_table = @embedFile("allocation_site_floors");

test "no family silently stops reaching allocation sites it used to reach" {
    for (families) |family| {
        errdefer std.debug.print("family: {s}\n", .{family.name});
        const measured = try core.oom.countAllocationSites(
            std.testing.allocator,
            generateAndDiscard,
            .{family},
        );
        try core.oom.expectSiteFloor(allocation_site_table, "api.generate", family.name, measured);
    }
}
