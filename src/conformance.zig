//! Conformance layer root. Oracle adapters, probes, and comparisons live
//! here; no production layer may import this module.

pub const probe_types = @import("conformance/probe_types.zig");
pub const mae = @import("conformance/mae.zig");
pub const corpus = @import("conformance/corpus.zig");
pub const comparison = @import("conformance/comparison.zig");

pub const DofProbe = probe_types.DofProbe;
pub const RingProbe = probe_types.RingProbe;
pub const TemplateMappingProbe = probe_types.TemplateMappingProbe;
pub const FragmentProbe = probe_types.FragmentProbe;
pub const ComponentProbe = probe_types.ComponentProbe;
pub const ComponentTransformStatus = probe_types.ComponentTransformStatus;
pub const dofProbe = probe_types.dofProbe;

test {
    _ = probe_types;
    _ = mae;
    _ = corpus;
    _ = comparison;
}
