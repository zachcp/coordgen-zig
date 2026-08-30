//! Seed tables for every fuzz target, gathered so a module can reach them by
//! one import name regardless of which layer it lives in.
//!
//! Each table is written by tools/run-fuzz, which proves a seed reproduces
//! before adding it. An empty table means nothing has been promoted yet - not
//! that promotion is unimplemented.
pub const c_abi_contract_seeds = @import("c_abi_contract_for_atom_and_bond_inputs.zig").seeds;
pub const c_abi_extended_input_seeds = @import("c_abi_extended_input_contract.zig").seeds;
pub const c_abi_ownership_sequence_seeds = @import("c_abi_ownership_sequences.zig").seeds;
