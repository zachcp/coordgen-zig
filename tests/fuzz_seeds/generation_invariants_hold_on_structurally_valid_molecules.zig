//! Promoted regression seeds for the `generationInvariants` target.
//!
//! Written by tools/run-fuzz, which proves a seed reproduces its failure
//! before adding it. Each entry replays under ordinary `zig build test`.
//!
//! Empty means no failure has been found and promoted yet - not that
//! promotion is unimplemented. `tools/run-fuzz --self-test` proves the
//! promotion path works against a deliberately planted failure.
pub const seeds: []const []const u8 = &.{};
