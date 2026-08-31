//! Bootstrap and promoted regression seeds for the `c abi contract for atom and bond inputs` target.
//!
//! The first entry is a deterministic valid C-input bootstrap. Subsequent
//! entries are promoted by tools/run-fuzz only after replay proves they
//! reproduce their failure under ordinary `zig build test`.
pub const seeds: []const []const u8 = &.{@import("c_abi_bootstrap.zig").valid_sequence};
