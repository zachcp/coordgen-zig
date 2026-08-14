const std = @import("std");
const core = @import("core");
const model = @import("model");
const topology = @import("topology");
const layout = @import("layout");
const optimize = @import("optimize");

pub const components = @import("components.zig");
pub const residues = @import("residues.zig");

test {
    _ = components;
    _ = residues;
}

/// Restricted proof result. Coordinates and both maps are allocator-owned;
/// coordinates are always serialized back to caller input order.
pub const Result = struct {
    allocator: std.mem.Allocator,
    coordinates: []core.math.Vec2,
    input_to_internal: []u32,
    internal_to_input: []core.ids.InputIndex,
    clean_pose: bool,

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
    if (prepared.working.active_atom_count != input.atoms.len) {
        return error.Unsupported;
    }
    var residue_atoms: []bool = &.{};
    defer if (residue_atoms.len != 0) allocator.free(residue_atoms);
    if (prepared.working.residues.len != 0) {
        residue_atoms = allocator.alloc(bool, prepared.working.atoms.len) catch return error.OutOfMemory;
        @memset(residue_atoms, false);
        for (prepared.working.residues) |residue| {
            if (residue.atom.index() >= residue_atoms.len or residue_atoms[residue.atom.index()] or
                prepared.graph.degree(residue.atom) != 0) return error.InvalidMapping;
            residue_atoms[residue.atom.index()] = true;
        }
    }
    const proximity_relations = try components.collectProximityRelations(allocator, prepared.working);
    defer allocator.free(proximity_relations);

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
    const clean_pose = try optimizeDiscrete(
        allocator,
        prepared.working.atoms,
        prepared.working.bonds,
        prepared.graph,
        prepared.rings,
        fragmentation,
        input.options.precision,
        residue_atoms,
    );
    try components.orientAcyclicComponents(
        allocator,
        prepared.working.atoms,
        prepared.working.bonds,
        prepared.graph,
        prepared.rings,
    );
    if (proximity_relations.len == 0 and residue_atoms.len == 0) {
        try components.arrangeComponents(allocator, prepared.working.atoms, prepared.graph);
    } else if (proximity_relations.len == 0) {
        try components.arrangeComponentsExcluding(allocator, prepared.working.atoms, prepared.graph, residue_atoms);
    } else if (prepared.working.residues.len != 0) {
        if (!try components.arrangeProximityComponentsExcluding(
            allocator,
            prepared.working.atoms,
            prepared.working.bonds,
            proximity_relations,
            residue_atoms,
        )) return error.Unsupported;
    } else if (!try components.arrangeProximityComponents(
        allocator,
        prepared.working.atoms,
        prepared.graph,
        prepared.rings,
        proximity_relations,
    )) {
        return error.Unsupported;
    }
    if (prepared.working.residues.len != 0) {
        if (prepared.working.residues.len == prepared.working.atoms.len) {
            try residues.placeProteinOnly(
                allocator,
                prepared.working.atoms,
                prepared.working.residues,
                prepared.working.string_bytes,
                prepared.working.residue_interactions,
            );
        } else {
            try residues.placeAroundLigand(
                allocator,
                prepared.working.atoms,
                prepared.working.bonds,
                prepared.working.residues,
                prepared.working.string_bytes,
                prepared.working.residue_interactions,
            );
        }
    }

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
        .clean_pose = clean_pose,
    };
}

const DiscreteScoreContext = struct {
    interactions: []const core.interaction.Interaction,
    bonds: []const optimize.discrete.BondScoreView,
    rings: []const optimize.discrete.RingScoreView,
    atom_fragments: []const core.ids.FragmentId,
    fragmentation: layout.Fragmentation,
};

fn scoreDiscretePose(raw_context: ?*anyopaque, coordinates: []const core.math.Vec2, dofs: []const core.dof.Dof) core.errors.Error!f32 {
    const context: *DiscreteScoreContext = @ptrCast(@alignCast(raw_context.?));
    var energy = try optimize.discrete.scoreClashInteractions(context.interactions, coordinates);
    energy += try optimize.discrete.scoreCrossBonds(context.bonds, coordinates);
    energy += try optimize.discrete.scoreAtomsInsideRings(context.rings, context.atom_fragments, coordinates);
    for (dofs) |dof| {
        const fragment = context.fragmentation.fragments[dof.fragment.index()];
        const chain_parent = fragment.parent.isValid() and fragment.flags.chain and
            context.fragmentation.fragments[fragment.parent.index()].flags.chain;
        energy += optimize.discrete.penalty(dof, dof.affected_atoms.len, .{
            .constrained_flip = fragment.flags.constrained_flip,
            .chain_with_chain_parent = chain_parent,
        }) catch |err| switch (err) {
            error.InvalidDofState => return error.InvalidMapping,
        };
    }
    return energy;
}

fn optimizeDiscrete(
    allocator: std.mem.Allocator,
    atoms: []model.Atom,
    bonds: []const model.Bond,
    graph: topology.Graph,
    rings: topology.RingMembership,
    fragmentation: layout.Fragmentation,
    precision: f32,
    excluded_atoms: []const bool,
) core.errors.Error!bool {
    if (excluded_atoms.len != 0 and excluded_atoms.len != atoms.len) return error.InvalidMapping;
    var dofs = try layout.macrocycle.collectAllDofs(allocator, bonds, graph, rings, fragmentation);
    defer dofs.deinit();
    if (dofs.items.len == 0) return false;

    const atom_has_dofs = allocator.alloc(bool, atoms.len) catch return error.OutOfMemory;
    defer allocator.free(atom_has_dofs);
    @memset(atom_has_dofs, false);
    for (dofs.affected_atoms) |atom| atom_has_dofs[atom.index()] = true;
    const atom_fragments = allocator.alloc(u32, atoms.len) catch return error.OutOfMemory;
    defer allocator.free(atom_fragments);
    for (fragmentation.atom_fragment, atom_fragments) |fragment, *output| output.* = fragment.index();
    var interactions = try optimize.buildBaseInteractions(allocator, atoms, bonds, .{
        .intrafragment_clashes = false,
        .atom_fragments = atom_fragments,
        .atom_has_dofs = atom_has_dofs,
    });
    defer interactions.deinit();
    var filtered_interactions: []core.interaction.Interaction = &.{};
    defer if (filtered_interactions.len != 0) allocator.free(filtered_interactions);
    var scoring_interactions = interactions.items;
    var filtered_count: usize = 0;
    if (excluded_atoms.len != 0) {
        filtered_interactions = allocator.alloc(core.interaction.Interaction, interactions.items.len) catch return error.OutOfMemory;
        for (interactions.items) |interaction| {
            if (interactionTouchesExcluded(interaction, excluded_atoms)) continue;
            filtered_interactions[filtered_count] = interaction;
            filtered_count += 1;
        }
        scoring_interactions = filtered_interactions[0..filtered_count];
    }

    const bond_views = allocator.alloc(optimize.discrete.BondScoreView, bonds.len) catch return error.OutOfMemory;
    defer allocator.free(bond_views);
    for (bonds, bond_views) |bond, *view| {
        var macrocycle = false;
        var small_ring = false;
        for (rings.bondRings(bond.id)) |ring_id| {
            const ring_size = rings.atoms(ring_id).len;
            macrocycle = macrocycle or ring_size >= topology.rings.macrocycle_size;
            small_ring = small_ring or ring_size < topology.rings.macrocycle_size;
        }
        view.* = .{
            .start = bond.start,
            .end = bond.end,
            .crossing_penalty_multiplier = bond.crossing_penalty_multiplier,
            .terminal = graph.degree(bond.start) == 1 or graph.degree(bond.end) == 1,
            .macrocycle = macrocycle,
            .small_ring = small_ring,
        };
    }
    const ring_views = allocator.alloc(optimize.discrete.RingScoreView, rings.rings.len) catch return error.OutOfMemory;
    defer allocator.free(ring_views);
    for (rings.rings, ring_views) |ring, *view| view.* = .{
        .atoms = rings.atoms(ring.id),
        .fragment = fragmentation.atom_fragment[rings.atoms(ring.id)[0].index()],
    };

    var frame_data = try layout.captureFragmentFrames(allocator, atoms, bonds, fragmentation);
    defer frame_data.deinit();
    const local_atoms = allocator.dupe(core.math.Vec2, frame_data.atom_coordinates) catch return error.OutOfMemory;
    defer allocator.free(local_atoms);
    const local_attachments = allocator.dupe(core.math.Vec2, frame_data.attachment_coordinates) catch return error.OutOfMemory;
    defer allocator.free(local_attachments);
    const global = allocator.alloc(core.math.Vec2, atoms.len) catch return error.OutOfMemory;
    defer allocator.free(global);
    var score_context = DiscreteScoreContext{
        .interactions = scoring_interactions,
        .bonds = bond_views,
        .rings = ring_views,
        .atom_fragments = fragmentation.atom_fragment,
        .fragmentation = fragmentation,
    };
    var evaluator = optimize.discrete.FramePoseEvaluator{
        .dofs = dofs.items,
        .affected_atoms = dofs.affected_atoms,
        .pose = .{
            .frames = frame_data.frames,
            .rebuild_order = frame_data.rebuild_order,
            .fragment_atoms = frame_data.fragment_atoms,
            .child_attachments = frame_data.child_attachments,
            .atom_coordinates = local_atoms,
            .attachment_coordinates = local_attachments,
            .global_coordinates = global,
        },
        .base_atom_coordinates = frame_data.atom_coordinates,
        .base_attachment_coordinates = frame_data.attachment_coordinates,
        .score_context = &score_context,
        .scorePoseFn = scoreDiscretePose,
    };
    var cache = try optimize.discrete.SolutionCache.init(allocator, dofs.items.len, precision, evaluator.evaluator());
    defer cache.deinit();
    const search = try optimize.discrete.tieredSearch(allocator, dofs.items, &cache, precision, 0);
    for (atoms, global) |*atom, coordinate| atom.coordinates = coordinate;
    return search.clean_pose;
}

fn interactionTouchesExcluded(interaction: core.interaction.Interaction, excluded_atoms: []const bool) bool {
    if (excluded_atoms.len == 0) return false;
    return switch (interaction.payload) {
        .stretch => |value| excluded_atoms[value.atom_a.index()] or excluded_atoms[value.atom_b.index()],
        .bend => |value| excluded_atoms[value.atom_a.index()] or excluded_atoms[value.center.index()] or excluded_atoms[value.atom_b.index()],
        .clash => |value| excluded_atoms[value.segment_start.index()] or excluded_atoms[value.point.index()] or excluded_atoms[value.segment_end.index()],
        .constraint => |value| excluded_atoms[value.atom.index()],
        .ez_constraint => |value| excluded_atoms[value.side_a.index()] or excluded_atoms[value.double_a.index()] or
            excluded_atoms[value.double_b.index()] or excluded_atoms[value.side_b.index()],
    };
}

fn rejectOutOfScope(input: anytype) core.errors.Error!void {
    if (input.extra_bonds.len != 0) return error.Unsupported;
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
