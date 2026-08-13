const std = @import("std");
const core = @import("core");
const topology = @import("topology");
const layout = @import("layout");

/// Restricted proof result. Coordinates and both maps are allocator-owned;
/// coordinates are always serialized back to caller input order.
pub const Result = struct {
    allocator: std.mem.Allocator,
    coordinates: []core.math.Vec2,
    input_to_internal: []u32,
    internal_to_input: []core.ids.InputIndex,

    pub fn deinit(self: *Result) void {
        self.allocator.free(self.internal_to_input);
        self.allocator.free(self.input_to_internal);
        self.allocator.free(self.coordinates);
        self.* = undefined;
    }
};

/// Compose validated safe-API-shaped input through native preparation,
/// fragmentation, and basic/macrocycle layout. This intentionally rejects
/// every domain whose owning native phase is not integrated yet.
pub fn generateValidated(allocator: std.mem.Allocator, input: anytype) core.errors.Error!Result {
    try rejectOutOfScope(input);
    var prepared = try topology.prepareInput(allocator, input);
    defer prepared.deinit();
    if (prepared.working.active_atom_count != input.atoms.len or
        prepared.graph.component_count != 1 or
        prepared.working.proximity_relations.len != 0)
    {
        return error.Unsupported;
    }

    var fragmentation = try layout.Fragmentation.init(
        allocator,
        prepared.working.atoms,
        prepared.working.bonds,
        prepared.graph,
        prepared.rings,
    );
    defer fragmentation.deinit();
    try layout.initializeCoordinatesWithOptions(
        allocator,
        prepared.working.atoms,
        prepared.working.bonds,
        prepared.graph,
        prepared.rings,
        fragmentation,
        input.options.force_open_macrocycles,
    );

    const coordinates = allocator.alloc(core.math.Vec2, input.atoms.len) catch return error.OutOfMemory;
    errdefer allocator.free(coordinates);
    const input_to_internal = allocator.dupe(u32, prepared.working.order.input_to_internal) catch return error.OutOfMemory;
    errdefer allocator.free(input_to_internal);
    const internal_to_input = allocator.dupe(core.ids.InputIndex, prepared.working.order.internal_to_input) catch return error.OutOfMemory;
    errdefer allocator.free(internal_to_input);
    for (prepared.working.atoms, prepared.working.order.internal_to_input) |atom, input_index| {
        coordinates[input_index] = atom.coordinates;
    }
    return .{
        .allocator = allocator,
        .coordinates = coordinates,
        .input_to_internal = input_to_internal,
        .internal_to_input = internal_to_input,
    };
}

fn rejectOutOfScope(input: anytype) core.errors.Error!void {
    if (input.residues.len != 0 or input.residue_interactions.len != 0 or input.extra_bonds.len != 0) return error.Unsupported;
    if (input.options.even_angles or input.options.constrain_all_atoms or input.options.build_from_fragments or
        input.options.debug_coordinates or input.options.template_directory != null) return error.Unsupported;
    for (input.atoms) |atom| {
        if (atom.hidden or atom.fixed or atom.constrained or atom.template_coordinates != null or
            atom.coordinates_3d != null or atom.stereo.value != .unspecified) return error.Unsupported;
    }
    for (input.bonds) |bond| {
        if (bond.skip or bond.stereo.value != .unspecified or bond.display != .none) return error.Unsupported;
    }
}
