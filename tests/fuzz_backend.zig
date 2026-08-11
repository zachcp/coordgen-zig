//! Refuses to let a fuzz gate be green while fuzzing nothing.
//!
//! The pinned toolchain's own fuzzer is the mechanism cgz-r15 chose, and
//! `std.testing.fuzz` in lib/compiler/test_runner.zig opens with:
//!
//!     // Some compiler backends are not capable of handling fuzz testing yet
//!     // but we still want CI test coverage enabled.
//!     if (need_simple) return;
//!
//! `need_simple` is true for a fixed set of self-hosted backends. On one of
//! those, every fuzz target returns immediately and the test binary reports
//! success - a gate that certifies nothing, which is worse than no gate. The
//! pin currently selects `stage2_llvm` on both CI targets, so fuzzing really
//! runs today; self-hosted aarch64 is under active development, so that is a
//! property of this toolchain version and not a safe assumption.
//!
//! The backend is read from `builtin.zig_backend` rather than inferred from
//! the target triple: which backend is default for a triple is exactly the
//! thing expected to change.
//!
//! The list below is a mirror of a private declaration, so it can go stale in
//! a direction that fails open - upstream adding a fourth non-fuzzing backend
//! would leave this list silently incomplete. tools/check-fuzz-backend reads
//! the switch out of the pinned toolchain's own test_runner.zig and fails if
//! the two disagree. See cgz-r27.

const std = @import("std");
const builtin = @import("builtin");

pub const Backend = @TypeOf(builtin.zig_backend);

/// Mirror of `need_simple` in the pinned toolchain's
/// lib/compiler/test_runner.zig. Keep sorted; tools/check-fuzz-backend
/// compares this set against that switch.
pub const non_fuzzing_backends = [_]Backend{
    .stage2_aarch64,
    .stage2_powerpc,
    .stage2_riscv64,
};

pub fn backendRunsFuzzTargets(backend: Backend) bool {
    for (non_fuzzing_backends) |non_fuzzing| {
        if (backend == non_fuzzing) return false;
    }
    return true;
}

/// Called by the fuzz harness before it claims to have fuzzed anything.
pub fn assertActiveBackendRunsFuzzTargets() error{FuzzTargetsAreNoOps}!void {
    if (!backendRunsFuzzTargets(builtin.zig_backend)) return error.FuzzTargetsAreNoOps;
}

test "the active compiler backend actually runs fuzz targets" {
    assertActiveBackendRunsFuzzTargets() catch {
        std.debug.print(
            \\
            \\The active compiler backend is {s}, for which the toolchain's
            \\std.testing.fuzz returns immediately. Every fuzz target would be
            \\a no-op reporting success. Select a backend that fuzzes (-fllvm
            \\on the pinned toolchain) or update cgz-r15's mechanism decision.
            \\
        , .{@tagName(builtin.zig_backend)});
        return error.FuzzTargetsAreNoOps;
    };
}

test "the predicate rejects exactly the mirrored backends" {
    // The host is one value; the predicate has to be right for all of them,
    // so it is exercised against the table rather than against this machine.
    for (non_fuzzing_backends) |backend| {
        try std.testing.expect(!backendRunsFuzzTargets(backend));
    }
    try std.testing.expect(backendRunsFuzzTargets(.stage2_llvm));
    try std.testing.expect(backendRunsFuzzTargets(.stage2_x86_64));
}
