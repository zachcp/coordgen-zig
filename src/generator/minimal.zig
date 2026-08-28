const std = @import("std");
const core = @import("core");
const model = @import("model");
const topology = @import("topology");
const layout = @import("layout");
const optimize = @import("optimize");

pub const components = @import("components.zig");
pub const residues = @import("residues.zig");
pub const bends = @import("bends.zig");
pub const inversions = @import("inversions.zig");

test {
    _ = components;
    _ = residues;
    _ = bends;
    _ = inversions;
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

/// Caller-owned destinations for one generation. Every slice is indexed by
/// caller input order; `coordinates`, `input_to_internal`, and `atom_stereo`
/// are atom-indexed, the two bond slices are bond-indexed, and
/// `internal_to_input` is the permutation itself. The three optional slices
/// exist because `generateValidated` does not need them; the public API
/// (cgz-7v2.21) does, and both go through one pipeline rather than two.
pub const Outputs = struct {
    coordinates: []core.math.Vec2,
    input_to_internal: []u32,
    internal_to_input: []core.ids.InputIndex,
    effective_bond_orders: ?[]core.chemistry.BondOrder = null,
    bond_displays: ?[]core.chemistry.BondDisplay = null,
    atom_stereo: ?[]core.chemistry.AtomStereo = null,
    /// The pose as it stands immediately before global orientation, in caller
    /// input order. Conformance-only: the oracle publishes the same stage
    /// through its probe hook, and comparing there separates a layout
    /// divergence from an orientation one, which the final coordinates alone
    /// cannot (cgz-7v2.4.2.1). No production caller requests it.
    pre_orientation: ?[]core.math.Vec2 = null,
    /// Conformance-only stage disposition. This makes the intentional
    /// fixed/constrained decline distinguishable from a silent input-class
    /// skip while leaving the public result unchanged (cgz-7v2.22).
    orientation_outcome: ?*components.OrientationOutcome = null,
    /// One generation-lifetime callback for the non-installed native probe.
    /// All values are borrowed and valid only for the callback invocation.
    probe_sink: ?ProbeSink = null,
};

pub const ProbeSnapshot = struct {
    prepared: *const topology.PreparedGraph,
    fragmentation: layout.Fragmentation,
    dofs: core.dof.Collection,
    clean_pose: bool,
};

pub const ProbeSink = struct {
    context: *anyopaque,
    template_observer: layout.TemplateObserver,
    captureFn: *const fn (*anyopaque, ProbeSnapshot) core.errors.Error!void,

    pub fn capture(self: ProbeSink, snapshot: ProbeSnapshot) core.errors.Error!void {
        return self.captureFn(self.context, snapshot);
    }
};

/// Compose validated safe-API-shaped input through native preparation,
/// fragmentation, and basic/macrocycle layout, writing every observable into
/// caller-owned storage and returning the clean-pose verdict. This
/// intentionally rejects every domain whose owning native phase is not
/// integrated yet.
///
/// Output storage is borrowed, never retained: on any error the destinations
/// hold unspecified values and the caller frees them, which is what lets the
/// public entry point own its result allocation in one place.
pub fn generateInto(allocator: std.mem.Allocator, input: anytype, outputs: Outputs) core.errors.Error!bool {
    if (outputs.coordinates.len != input.atoms.len or
        outputs.input_to_internal.len != input.atoms.len or
        outputs.internal_to_input.len != input.atoms.len)
    {
        return error.InvalidMapping;
    }
    if (outputs.atom_stereo) |slice| {
        if (slice.len != input.atoms.len) return error.InvalidMapping;
    }
    if (outputs.effective_bond_orders) |slice| {
        if (slice.len != input.bonds.len) return error.InvalidMapping;
    }
    if (outputs.bond_displays) |slice| {
        if (slice.len != input.bonds.len) return error.InvalidMapping;
    }
    if (outputs.pre_orientation) |slice| {
        if (slice.len != input.atoms.len) return error.InvalidMapping;
    }
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
    // One collection per generation, shared by layout and the discrete pass.
    // `collectAllDofs` reads only bonds, graph, ring membership and
    // fragmentation - no coordinates - so collecting it before layout is
    // valid, and layout needs it for the stages upstream gates on
    // `getDofsOfAtom(a).empty()` (cgz-7v2.28).
    var dofs = try layout.macrocycle.collectAllDofs(
        allocator,
        prepared.working.atoms,
        prepared.working.bonds,
        prepared.graph,
        prepared.rings,
        fragmentation,
    );
    defer dofs.deinit();
    const layout_outcome = try layout.initializeCoordinatesWithOptions(
        allocator,
        prepared.working.atoms,
        prepared.working.bonds,
        prepared.graph,
        prepared.rings,
        fragmentation,
        dofs,
        .{
            .force_open_macrocycles = input.options.force_open_macrocycles,
            .templates = .{
                .load_templates = input.options.load_templates,
                .template_directory = input.options.template_directory,
            },
            .template_observer = if (outputs.probe_sink) |sink| sink.template_observer else null,
        },
    );
    const clean_pose = try optimizeDiscrete(
        allocator,
        prepared.working.atoms,
        prepared.working.bonds,
        prepared.graph,
        prepared.rings,
        fragmentation,
        dofs,
        proximity_relations,
        input.options.precision,
        residue_atoms,
    );
    // Upstream requires minimization when the discrete search could not reach
    // a clean pose, and separately when maybeMinimizeRings fired during ring
    // placement (CoordgenMinimizer.cpp:1272 and :1458). Both feed one
    // molecule-level flag, and skip_minimization suppresses the run itself
    // rather than the flag, matching CoordgenMinimizer::run.
    if ((!clean_pose or layout_outcome.minimization_required) and
        !input.options.skip_minimization)
    {
        try minimizeGenerated(
            allocator,
            prepared.working.atoms,
            prepared.working.bonds,
            prepared.graph,
            prepared.rings,
            fragmentation,
            input.options.even_angles,
        );
    }
    if (outputs.pre_orientation) |slice| {
        for (prepared.working.atoms, prepared.working.order.internal_to_input) |atom, input_index| {
            if (input_index >= slice.len) return error.InvalidMapping;
            slice[input_index] = atom.coordinates;
        }
    }
    const orientation_outcome = try components.orientComponents(
        allocator,
        prepared.working.atoms,
        prepared.working.bonds,
        prepared.graph,
        prepared.rings,
        fragmentation,
    );
    if (outputs.orientation_outcome) |outcome| outcome.* = orientation_outcome;
    if (proximity_relations.len == 0 and residue_atoms.len == 0) {
        try components.arrangeComponents(allocator, prepared.working.atoms, prepared.graph);
    } else if (proximity_relations.len == 0) {
        try components.arrangeComponentsExcluding(allocator, prepared.working.atoms, prepared.graph, residue_atoms);
    } else {
        try components.arrangeComponentsWithProximityExcluding(
            allocator,
            prepared.working.atoms,
            prepared.working.bonds,
            prepared.graph,
            proximity_relations,
            residue_atoms,
        );
    }
    if (outputs.probe_sink) |sink| try sink.capture(.{
        .prepared = &prepared,
        .fragmentation = fragmentation,
        .dofs = dofs,
        .clean_pose = clean_pose,
    });
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
    for (prepared.working.atoms) |atom| if (!atom.coordinates.isFinite()) return error.InvalidCoordinate;
    try topology.stereo.writeAtomBondDisplays(
        allocator,
        prepared.working.atoms,
        prepared.working.bonds,
        prepared.graph,
        prepared.rings,
    );

    if (prepared.working.order.input_to_internal.len != outputs.input_to_internal.len or
        prepared.working.order.internal_to_input.len != outputs.internal_to_input.len)
    {
        return error.InvalidMapping;
    }
    @memcpy(outputs.input_to_internal, prepared.working.order.input_to_internal);
    @memcpy(outputs.internal_to_input, prepared.working.order.internal_to_input);
    for (prepared.working.atoms, prepared.working.order.internal_to_input) |atom, input_index| {
        outputs.coordinates[input_index] = atom.coordinates;
    }
    if (outputs.atom_stereo) |slice| {
        for (prepared.working.atoms) |atom| {
            if (atom.input_index >= slice.len) return error.InvalidMapping;
            slice[atom.input_index] = atom.stereo;
        }
    }
    // Bonds are stored in structural-first internal order, so they are
    // serialized through their own input index rather than the atom
    // permutation. prepare() emits exactly one internal bond per input bond.
    if (outputs.effective_bond_orders) |slice| {
        if (prepared.working.bonds.len != slice.len) return error.InvalidMapping;
        for (prepared.working.bonds) |bond| {
            if (bond.input_index >= slice.len) return error.InvalidMapping;
            slice[bond.input_index] = bond.effective_order;
        }
    }
    if (outputs.bond_displays) |slice| {
        if (prepared.working.bonds.len != slice.len) return error.InvalidMapping;
        for (prepared.working.bonds) |bond| {
            if (bond.input_index >= slice.len) return error.InvalidMapping;
            slice[bond.input_index] = bond.display;
        }
    }
    return clean_pose;
}

/// Restricted proof entry point: allocates and owns the three observables it
/// reports. Retained for the native module tests that predate the public
/// entry point; new callers should use the public `api.generate`.
pub fn generateValidated(allocator: std.mem.Allocator, input: anytype) core.errors.Error!Result {
    const coordinates = allocator.alloc(core.math.Vec2, input.atoms.len) catch return error.OutOfMemory;
    errdefer allocator.free(coordinates);
    const input_to_internal = allocator.alloc(u32, input.atoms.len) catch return error.OutOfMemory;
    errdefer allocator.free(input_to_internal);
    const internal_to_input = allocator.alloc(core.ids.InputIndex, input.atoms.len) catch return error.OutOfMemory;
    errdefer allocator.free(internal_to_input);
    const clean_pose = try generateInto(allocator, input, .{
        .coordinates = coordinates,
        .input_to_internal = input_to_internal,
        .internal_to_input = internal_to_input,
    });
    return .{
        .allocator = allocator,
        .coordinates = coordinates,
        .input_to_internal = input_to_internal,
        .internal_to_input = internal_to_input,
        .clean_pose = clean_pose,
    };
}

/// Assemble the four interaction families upstream feeds
/// addInteractionsOfMolecule - clash and stretch with intrafragment clashes
/// enabled, bends, and chiral inversion constraints - and run the continuous
/// minimizer over them.
///
/// The clash set differs from the discrete pass on purpose: minimizeMolecule
/// passes intrafragmentClashes = true where avoidClashesOfMolecule passes
/// false.
fn minimizeGenerated(
    allocator: std.mem.Allocator,
    atoms: []model.Atom,
    bonds: []const model.Bond,
    graph: topology.Graph,
    rings: topology.RingMembership,
    fragmentation: layout.Fragmentation,
    even_angles: bool,
) core.errors.Error!void {
    const atom_fragments = allocator.alloc(u32, atoms.len) catch return error.OutOfMemory;
    defer allocator.free(atom_fragments);
    for (fragmentation.atom_fragment, atom_fragments) |fragment, *output| output.* = fragment.index();

    var base = try optimize.buildBaseInteractions(allocator, atoms, bonds, .{
        .intrafragment_clashes = true,
        .atom_fragments = atom_fragments,
    });
    defer base.deinit();

    var analysis = try topology.rings.Analysis.init(allocator, rings, atoms, bonds);
    defer analysis.deinit();
    var groups = try bends.build(allocator, atoms, bonds, graph, rings, analysis);
    defer groups.deinit();
    var bend_interactions = try optimize.buildBendInteractions(
        allocator,
        atoms,
        groups.groups,
        .{ .even_angles = even_angles },
    );
    defer bend_interactions.deinit();

    const ez_constraints = try inversions.buildChiralInversionConstraints(
        allocator,
        atoms,
        bonds,
        graph,
        rings,
    );
    defer allocator.free(ez_constraints);
    var ez_interactions = try optimize.buildEzInteractions(allocator, ez_constraints);
    defer ez_interactions.deinit();

    var combined = try optimize.combineInteractions(allocator, &.{
        base.items,
        bend_interactions.items,
        ez_interactions.items,
    });
    defer combined.deinit();

    _ = try optimize.minimizeMolecule(allocator, atoms, combined.items, .{}, null);
}

const DiscreteScoreContext = struct {
    interactions: []const core.interaction.Interaction,
    bonds: []const optimize.discrete.BondScoreView,
    rings: []const optimize.discrete.RingScoreView,
    atom_fragments: []const core.ids.FragmentId,
    fragment_components: []const core.ids.MoleculeId,
    fragmentation: layout.Fragmentation,
    /// Upstream's proximity term (below). Empty whenever the molecule has no
    /// proximity relations, which is the common case and costs nothing.
    proximity: Proximity = .{},

    const Proximity = struct {
        sites: []const components.ProximityScoreSite = &.{},
        /// Refilled per scored pose; `sites` supplies everything but the
        /// addition vector.
        views: []optimize.discrete.ProximityScoreView = &.{},
        /// The addition vector is defined on atoms, and the search scores bare
        /// coordinate sets, so this holds the candidate pose in the shape
        /// `components.singleAdditionVector` reads. It is a scratch copy: the
        /// alternative was a second coordinate-taking implementation of the
        /// same eight lines, which is the duplication cgz-7v2.24 exists to
        /// find.
        atoms: []model.Atom = &.{},
        graph: ?topology.Graph = null,
        rings: ?topology.RingMembership = null,
    };
};

fn scoreDiscretePose(raw_context: ?*anyopaque, coordinates: []const core.math.Vec2, dofs: []const core.dof.Dof) core.errors.Error!f32 {
    const context: *DiscreteScoreContext = @ptrCast(@alignCast(raw_context.?));
    var energy = try optimize.discrete.scoreClashInteractions(context.interactions, coordinates);
    energy += try optimize.discrete.scoreCrossBonds(context.bonds, coordinates);
    energy += try optimize.discrete.scoreAtomsInsideRings(context.rings, context.atom_fragments, context.fragment_components, coordinates);
    // `CoordgenMinimizer::scoreClashes` ends with this term
    // (CoordgenMinimizer.cpp:808). Its `scoreProximityRelationsOnOppositeSides`
    // parameter defaults true and no caller ever passes false - :1268 passes a
    // local named `doNotComputeForces` that is initialised to true - so it is
    // unconditional in practice. It penalises a pose in which two of a
    // molecule's proximity attachment points face more than 90 degrees apart,
    // and it is purely intramolecular: the addition vectors are offsets from
    // an atom to its own neighbours, so this does not depend on where the
    // molecules have been arranged relative to each other.
    if (context.proximity.sites.len != 0) {
        const proximity = context.proximity;
        if (proximity.atoms.len != coordinates.len) return error.InvalidMapping;
        for (proximity.atoms, coordinates) |*atom, coordinate| atom.coordinates = coordinate;
        for (proximity.sites, proximity.views) |site, *view| view.* = .{
            .local_molecule = site.local_molecule,
            .local_fragment = site.local_fragment,
            .other_molecule = site.other_molecule,
            .addition_vector = components.singleAdditionVector(
                proximity.atoms,
                proximity.graph.?,
                proximity.rings.?,
                site.local_atom,
            ),
        };
        energy += try optimize.discrete.scoreProximityRelationsOnOppositeSides(proximity.views);
    }
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
    dofs: core.dof.Collection,
    proximity_relations: []const components.ProximityRelation,
    precision: f32,
    excluded_atoms: []const bool,
) core.errors.Error!bool {
    if (excluded_atoms.len != 0 and excluded_atoms.len != atoms.len) return error.InvalidMapping;
    const atom_has_dofs = allocator.alloc(bool, atoms.len) catch return error.OutOfMemory;
    defer allocator.free(atom_has_dofs);
    @memset(atom_has_dofs, false);
    for (dofs.affected_atoms) |atom| atom_has_dofs[atom.index()] = true;
    const atom_fragments = allocator.alloc(u32, atoms.len) catch return error.OutOfMemory;
    defer allocator.free(atom_fragments);
    for (fragmentation.atom_fragment, atom_fragments) |fragment, *output| output.* = fragment.index();
    var base = try optimize.buildBaseInteractions(allocator, atoms, bonds, .{
        .intrafragment_clashes = false,
        .atom_fragments = atom_fragments,
        .atom_has_dofs = atom_has_dofs,
    });
    defer base.deinit();
    // Upstream's discrete pass is addClashInteractionsOfMolecule followed by
    // addPeptideBondInversionConstraintsOfMolecule, and only then scoreClashes
    // and flipFragments (CoordgenMinimizer.cpp:1254-1263). The peptide
    // constraints therefore influence the clean-pose verdict, not just the
    // geometry.
    const peptide_constraints = try inversions.buildPeptideBondInversionConstraints(
        allocator,
        atoms,
        bonds,
        graph,
    );
    defer allocator.free(peptide_constraints);
    var peptide_interactions = try optimize.buildEzInteractions(allocator, peptide_constraints);
    defer peptide_interactions.deinit();
    var interactions = try optimize.combineInteractions(allocator, &.{
        base.items,
        peptide_interactions.items,
    });
    defer interactions.deinit();
    var filtered_interactions: []core.interaction.Interaction = &.{};
    defer if (filtered_interactions.len != 0) allocator.free(filtered_interactions);
    var scoring_interactions = interactions.items;
    var filtered_count: usize = 0;
    if (excluded_atoms.len != 0 or graph.component_count > 1) {
        filtered_interactions = allocator.alloc(core.interaction.Interaction, interactions.items.len) catch return error.OutOfMemory;
        for (interactions.items) |interaction| {
            if (interactionTouchesExcluded(interaction, excluded_atoms)) continue;
            if (try interactionCrossesComponents(interaction, graph)) continue;
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
            .component = graph.component(bond.start) orelse return error.InvalidMapping,
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
        .component = graph.component(rings.atoms(ring.id)[0]) orelse return error.InvalidMapping,
    };
    const fragment_components = allocator.alloc(core.ids.MoleculeId, fragmentation.fragments.len) catch return error.OutOfMemory;
    defer allocator.free(fragment_components);
    for (fragmentation.fragments, fragment_components) |fragment, *component| component.* = fragment.component;
    const proximity_sites = try components.proximityScoreSites(
        allocator,
        graph,
        fragmentation.atom_fragment,
        proximity_relations,
    );
    defer allocator.free(proximity_sites);
    var all_clean = true;
    for (0..graph.component_count) |raw_component| {
        const component = core.ids.MoleculeId.fromIndex(@intCast(raw_component));
        const clean = try optimizeDiscreteComponent(
            allocator,
            atoms,
            graph,
            rings,
            fragmentation,
            dofs,
            scoring_interactions,
            bond_views,
            ring_views,
            fragment_components,
            proximity_sites,
            component,
            precision,
        );
        all_clean = all_clean and clean;
    }
    return all_clean;
}

fn optimizeDiscreteComponent(
    allocator: std.mem.Allocator,
    atoms: []model.Atom,
    graph: topology.Graph,
    rings: topology.RingMembership,
    fragmentation: layout.Fragmentation,
    dofs: core.dof.Collection,
    all_interactions: []const core.interaction.Interaction,
    all_bonds: []const optimize.discrete.BondScoreView,
    all_rings: []const optimize.discrete.RingScoreView,
    fragment_components: []const core.ids.MoleculeId,
    all_proximity_sites: []const components.ProximityScoreSite,
    component: core.ids.MoleculeId,
    precision: f32,
) core.errors.Error!bool {
    var component_dofs: std.ArrayList(core.dof.Dof) = .empty;
    defer component_dofs.deinit(allocator);
    var dof_indices: std.ArrayList(usize) = .empty;
    defer dof_indices.deinit(allocator);
    for (dofs.items, 0..) |dof, index| {
        if (fragmentation.fragments[dof.fragment.index()].component != component) continue;
        component_dofs.append(allocator, dof) catch return error.OutOfMemory;
        dof_indices.append(allocator, index) catch return error.OutOfMemory;
    }

    var interactions: std.ArrayList(core.interaction.Interaction) = .empty;
    defer interactions.deinit(allocator);
    for (all_interactions) |interaction| {
        if (try interactionComponent(interaction, graph) == component) {
            interactions.append(allocator, interaction) catch return error.OutOfMemory;
        }
    }
    var bonds: std.ArrayList(optimize.discrete.BondScoreView) = .empty;
    defer bonds.deinit(allocator);
    for (all_bonds) |bond| if (bond.component == component) {
        bonds.append(allocator, bond) catch return error.OutOfMemory;
    };
    var proximity_sites: std.ArrayList(components.ProximityScoreSite) = .empty;
    defer proximity_sites.deinit(allocator);
    for (all_proximity_sites) |site| if (site.local_molecule == component) {
        proximity_sites.append(allocator, site) catch return error.OutOfMemory;
    };

    var score_context = DiscreteScoreContext{
        .interactions = interactions.items,
        .bonds = bonds.items,
        // Upstream's atom-inside-ring term is global even while the mutable
        // clash/DOF search is molecule-local. Each ring still compares only
        // with atoms from its own molecule in scoreAtomsInsideRings.
        .rings = all_rings,
        .atom_fragments = fragmentation.atom_fragment,
        .fragment_components = fragment_components,
        .fragmentation = fragmentation,
    };
    var proximity_views: []optimize.discrete.ProximityScoreView = &.{};
    defer if (proximity_views.len != 0) allocator.free(proximity_views);
    var proximity_atoms: []model.Atom = &.{};
    defer if (proximity_atoms.len != 0) allocator.free(proximity_atoms);
    if (proximity_sites.items.len != 0) {
        proximity_views = allocator.alloc(optimize.discrete.ProximityScoreView, proximity_sites.items.len) catch return error.OutOfMemory;
        proximity_atoms = allocator.dupe(model.Atom, atoms) catch return error.OutOfMemory;
        score_context.proximity = .{
            .sites = proximity_sites.items,
            .views = proximity_views,
            .atoms = proximity_atoms,
            .graph = graph,
            .rings = rings,
        };
    }

    const coordinates = allocator.alloc(core.math.Vec2, atoms.len) catch return error.OutOfMemory;
    defer allocator.free(coordinates);
    for (atoms, coordinates) |atom, *coordinate| coordinate.* = atom.coordinates;
    const initial = try scoreDiscretePose(&score_context, coordinates, component_dofs.items);
    if (initial < optimize.discrete.clash_energy_threshold) return true;
    if (component_dofs.items.len == 0) return false;

    var frame_data = try layout.captureFragmentFrames(allocator, atoms, fragmentation);
    defer frame_data.deinit();
    const local_atoms = allocator.dupe(core.math.Vec2, frame_data.atom_coordinates) catch return error.OutOfMemory;
    defer allocator.free(local_atoms);
    const local_attachments = allocator.dupe(core.math.Vec2, frame_data.attachment_coordinates) catch return error.OutOfMemory;
    defer allocator.free(local_attachments);
    const global = allocator.alloc(core.math.Vec2, atoms.len) catch return error.OutOfMemory;
    defer allocator.free(global);
    var evaluator = optimize.discrete.FramePoseEvaluator{
        .dofs = component_dofs.items,
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
    var cache = try optimize.discrete.SolutionCache.init(allocator, component_dofs.items.len, precision, evaluator.evaluator());
    defer cache.deinit();
    const search = try optimize.discrete.tieredSearch(allocator, component_dofs.items, &cache, precision, 0);
    for (component_dofs.items, dof_indices.items) |source, index| dofs.items[index].state = source.state;
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

fn interactionCrossesComponents(interaction: core.interaction.Interaction, graph: topology.Graph) core.errors.Error!bool {
    const pair = switch (interaction.payload) {
        .clash => |value| .{ value.segment_start, value.point },
        else => return false,
    };
    const first = graph.component(pair[0]) orelse return error.InvalidMapping;
    const second = graph.component(pair[1]) orelse return error.InvalidMapping;
    return first != second;
}

fn interactionComponent(interaction: core.interaction.Interaction, graph: topology.Graph) core.errors.Error!core.ids.MoleculeId {
    const atom = switch (interaction.payload) {
        .stretch => |value| value.atom_a,
        .bend => |value| value.center,
        .clash => |value| value.segment_start,
        .constraint => |value| value.atom,
        .ez_constraint => |value| value.double_a,
    };
    return graph.component(atom) orelse error.InvalidMapping;
}

fn rejectOutOfScope(input: anytype) core.errors.Error!void {
    if (input.extra_bonds.len != 0) return error.Unsupported;
    if (input.options.even_angles or input.options.constrain_all_atoms or
        input.options.debug_coordinates or input.options.template_directory != null) return error.Unsupported;
    for (input.atoms) |atom| {
        if (atom.hidden or atom.fixed or atom.constrained or atom.template_coordinates != null or
            atom.coordinates_3d != null or atom.stereo.value != .unspecified) return error.Unsupported;
    }
    for (input.bonds) |bond| {
        if (bond.skip or bond.stereo.value != .unspecified or bond.display != .none) return error.Unsupported;
    }
}
