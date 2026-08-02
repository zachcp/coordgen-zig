const std = @import("std");

/// Import direction is `from` imports `to`. The build graph must install this
/// exact allow-list with explicit addImport calls; Zig itself permits cycles.
pub const Layer = enum(u8) {
    core,
    model,
    geometry,
    topology,
    layout,
    optimize,
    generator,
    api,
    c_abi,
    conformance,
};

pub const Edge = struct {
    from: Layer,
    to: Layer,
};

pub const approved_edges = [_]Edge{
    .{ .from = .model, .to = .core },
    .{ .from = .geometry, .to = .core },
    .{ .from = .topology, .to = .core },
    .{ .from = .topology, .to = .model },
    .{ .from = .topology, .to = .geometry },
    .{ .from = .layout, .to = .core },
    .{ .from = .layout, .to = .model },
    .{ .from = .layout, .to = .geometry },
    .{ .from = .layout, .to = .topology },
    .{ .from = .optimize, .to = .core },
    .{ .from = .optimize, .to = .model },
    .{ .from = .optimize, .to = .geometry },
    .{ .from = .generator, .to = .core },
    .{ .from = .generator, .to = .model },
    .{ .from = .generator, .to = .geometry },
    .{ .from = .generator, .to = .topology },
    .{ .from = .generator, .to = .layout },
    .{ .from = .generator, .to = .optimize },
    .{ .from = .api, .to = .core },
    .{ .from = .api, .to = .generator },
    .{ .from = .c_abi, .to = .core },
    .{ .from = .c_abi, .to = .api },
    .{ .from = .conformance, .to = .core },
    .{ .from = .conformance, .to = .model },
    .{ .from = .conformance, .to = .geometry },
    .{ .from = .conformance, .to = .topology },
    .{ .from = .conformance, .to = .layout },
    .{ .from = .conformance, .to = .optimize },
    .{ .from = .conformance, .to = .generator },
    .{ .from = .conformance, .to = .api },
};

pub fn allows(from: Layer, to: Layer) bool {
    for (approved_edges) |edge| {
        if (edge.from == from and edge.to == to) return true;
    }
    return false;
}

fn rank(layer: Layer) u8 {
    return switch (layer) {
        .core => 0,
        .model, .geometry => 1,
        .topology => 2,
        .layout, .optimize => 3,
        .generator => 4,
        .api => 5,
        .c_abi => 6,
        .conformance => 7,
    };
}

pub fn validate() !void {
    for (approved_edges, 0..) |edge, index| {
        if (edge.from == edge.to or rank(edge.from) <= rank(edge.to)) {
            return error.InvalidModuleEdge;
        }
        if (edge.to == .conformance) return error.ProductionImportsConformance;
        for (approved_edges[index + 1 ..]) |other| {
            if (edge.from == other.from and edge.to == other.to) {
                return error.DuplicateModuleEdge;
            }
        }
    }
}

test "approved module graph is one-way and production cannot import conformance" {
    try validate();
    try std.testing.expect(allows(.model, .core));
    try std.testing.expect(allows(.optimize, .core));
    try std.testing.expect(!allows(.optimize, .layout));
    try std.testing.expect(!allows(.generator, .conformance));
}
