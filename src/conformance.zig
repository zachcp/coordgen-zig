//! Conformance layer root. Oracle adapters, probes, and comparisons live
//! here; no production layer may import this module.

pub const probe_types = @import("conformance/probe_types.zig");
pub const mae = @import("conformance/mae.zig");

pub const DofProbe = probe_types.DofProbe;
pub const dofProbe = probe_types.dofProbe;

test {
    _ = probe_types;
    _ = mae;
}
