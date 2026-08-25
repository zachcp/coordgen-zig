//! Seed tables for every fuzz target, gathered so a module can reach them by
//! one import name regardless of which layer it lives in.
//!
//! Each table is written by tools/run-fuzz, which proves a seed reproduces
//! before adding it. An empty table means nothing has been promoted yet - not
//! that promotion is unimplemented.
pub const c_abi_contract_seeds = @import("c_abi_contract_holds_for_any_input_a_c_caller_can_express.zig").seeds;
