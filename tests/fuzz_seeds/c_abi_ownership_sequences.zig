//! Bootstrap and promoted regression seeds for the `c abi ownership sequences` target.
//!
//! The first entry is a deterministic bootstrap corpus input that reaches a
//! valid multi-result ownership sequence even from a fresh fuzzer cache.
//! Subsequent entries are promoted by tools/run-fuzz only after replay proves
//! they reproduce their failure under ordinary `zig build test`.
pub const seeds: []const []const u8 = &.{@import("c_abi_bootstrap.zig").valid_sequence};
