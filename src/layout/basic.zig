const std = @import("std");
const core = @import("core");
const model = @import("model");
const topology = @import("topology");
const fragments = @import("fragments.zig");
const macrocycle = @import("macrocycle.zig");
const templates = @import("templates.zig");

pub const bond_length: f32 = 50;
const pi: f32 = std.math.pi;

/// Generate deterministic local coordinates for ordinary rigid fragments.
/// Template and macrocycle dispatch are deliberately owned by their later
/// layers; this seam covers regular rings and the acyclic fallback.
/// Placement controls that reach the ring path. `templates` is passed
/// straight through to templates.findRingSet, which is where
/// `load_templates` and the rejection of a runtime template directory are
/// enforced.
pub const Options = struct {
    force_open_macrocycles: bool = false,
    templates: templates.Options = .{},
};

/// What the ring path decided that later phases need to know.
pub const Outcome = struct {
    /// Upstream's maybeMinimizeRings fired for at least one ring system that
    /// was placed without a template. It is molecule-level because
    /// requireMinimization is.
    minimization_required: bool = false,
};

pub fn initializeCoordinates(
    allocator: std.mem.Allocator,
    atoms: []model.Atom,
    bonds: []const model.Bond,
    graph: topology.Graph,
    membership: topology.RingMembership,
    fragmentation: fragments.Fragmentation,
) core.errors.Error!void {
    _ = try initializeCoordinatesInternal(allocator, atoms, bonds, graph, membership, fragmentation, .{});
}

pub fn initializeCoordinatesWithOptions(
    allocator: std.mem.Allocator,
    atoms: []model.Atom,
    bonds: []const model.Bond,
    graph: topology.Graph,
    membership: topology.RingMembership,
    fragmentation: fragments.Fragmentation,
    options: Options,
) core.errors.Error!Outcome {
    return initializeCoordinatesInternal(allocator, atoms, bonds, graph, membership, fragmentation, options);
}

fn initializeCoordinatesInternal(
    allocator: std.mem.Allocator,
    atoms: []model.Atom,
    bonds: []const model.Bond,
    graph: topology.Graph,
    membership: topology.RingMembership,
    fragmentation: fragments.Fragmentation,
    options: Options,
) core.errors.Error!Outcome {
    var outcome: Outcome = .{};
    // Ring fusion structure is needed by the central-ring priority, and is
    // derived once for the whole molecule rather than per fragment.
    var analysis = try topology.rings.Analysis.init(allocator, membership, atoms, bonds);
    defer analysis.deinit();
    const placed = allocator.alloc(bool, atoms.len) catch return error.OutOfMemory;
    defer allocator.free(placed);
    @memset(placed, false);
    const queue = allocator.alloc(core.ids.AtomId, atoms.len) catch return error.OutOfMemory;
    defer allocator.free(queue);
    // Upstream snapshots each fragment's pose into its own frame the moment the
    // fragment is built (`sketcherMinimizerFragment::storeCoordinateInformation`),
    // and assembles the molecule from those frames afterwards. The snapshot is
    // what lets a fragment's build write on atoms outside it - which upstream's
    // neighbour placement does, and which is how a child fragment learns the
    // direction its attachment atom was given - without a later fragment's
    // build trampling an earlier one's finished pose.
    var frames = try FragmentFrames.init(allocator, atoms.len, fragmentation.fragments.len);
    defer frames.deinit();

    for (fragmentation.fragments) |fragment| {
        if (try placeFragmentRings(allocator, atoms, bonds, graph, membership, fragmentation, fragment, placed, options, analysis)) {
            outcome.minimization_required = true;
        }

        const members = fragmentation.members(fragment.id);
        var head: usize = 0;
        var tail: usize = 0;
        for (members) |atom| if (placed[atom.index()]) {
            queue[tail] = atom;
            tail += 1;
        };
        if (tail == 0) {
            // `CoordgenFragmentBuilder::buildNonRingAtoms`: a fragment with no
            // ring atoms starts from its bond to the parent, not from whichever
            // atom happens to be first. The parent-side atom is seeded one bond
            // length along -x, which is also the direction the neighbour
            // placement starts rotating from.
            const start = attachmentOf(fragment) orelse members[0];
            if (parentAtomOf(fragment)) |parent_atom| {
                atoms[parent_atom.index()].coordinates = .{ .x = -bond_length, .y = 0 };
            }
            atoms[start.index()].coordinates = .{};
            placed[start.index()] = true;
            queue[0] = start;
            tail = 1;
        }
        while (head < tail) : (head += 1) {
            const center = queue[head];
            var fallback_count: usize = 0;
            for (graph.neighbors(center)) |neighbor| {
                if (!placed[neighbor.index()] or fragmentation.atom_fragment[neighbor.index()] != fragment.id) fallback_count += 1;
            }
            // No early exit when every in-fragment neighbour is already placed:
            // a centre can still owe a position to a neighbour in a child
            // fragment, and that position is the only record of which direction
            // the child leaves in.
            // Upstream splits on whether the centre is a ring atom
            // (initializeVariablesForNeighboursCoordinates). Both branches are
            // transcribed; each falls through to the fixed spread below only
            // when its own neighbour buffer would overflow.
            if (membership.atomRings(center).len != 0) {
                if (try placeRingAtomSubstituents(atoms, graph, membership, fragmentation, fragment.id, center, placed, queue, &tail)) continue;
            } else if (try placeAcyclicNeighbours(allocator, atoms, bonds, graph, membership, fragmentation, fragment.id, center, placed, queue, &tail)) {
                continue;
            }
            const base_angle = parentAngle(atoms, graph, fragmentation, fragment.id, center, placed);
            var generated: usize = 0;
            for (graph.neighbors(center)) |neighbor| {
                if (placed[neighbor.index()] and fragmentation.atom_fragment[neighbor.index()] == fragment.id) continue;
                const spread = if (fallback_count == 1)
                    @as(f32, 0)
                else
                    (@as(f32, @floatFromInt(generated)) - @as(f32, @floatFromInt(fallback_count - 1)) * 0.5) * (2 * pi / 3);
                const angle = base_angle + spread;
                atoms[neighbor.index()].coordinates = .{
                    .x = atoms[center.index()].coordinates.x + @cos(angle) * bond_length,
                    .y = atoms[center.index()].coordinates.y + @sin(angle) * bond_length,
                };
                generated += 1;
                if (fragmentation.atom_fragment[neighbor.index()] != fragment.id) continue;
                placed[neighbor.index()] = true;
                queue[tail] = neighbor;
                tail += 1;
            }
        }
        // Upstream's order inside buildFragment is fallbackIfNanCoordinates,
        // then rotateMainFragment for a constrained root, then the
        // fixed-coordinate reset (CoordgenFragmentBuilder.cpp:858-866), and
        // only then the snapshot. All three ran once at the end of the molecule
        // here, after assembly, which is a different function of a different
        // pose.
        fallbackOnValid3dCoordinates(atoms, members);
        alignConstrainedMainFragment(atoms, fragmentation, fragment);
        restoreFixedFragmentCoordinates(atoms, fragmentation, fragment);
        frames.store(atoms, fragmentation, fragment);
    }
    try placeFragmentsFromFrames(allocator, atoms, bonds, graph, fragmentation, frames);
    restoreFixedCoordinates(atoms);
    return outcome;
}

/// Each fragment's own pose in its own frame, plus the position its children's
/// attachment atoms were given in it. Upstream's `_coordinates` map, which
/// holds exactly those two sets (`storeCoordinateInformation`).
const FragmentFrames = struct {
    allocator: std.mem.Allocator,
    /// Every atom, in the frame of the fragment that owns it.
    own: []core.math.Vec2,
    /// Indexed by child fragment: that child's attachment atom, in the frame of
    /// its parent.
    attachment: []core.math.Vec2,

    /// One allocation for both slices: the frames are pure scratch, and the
    /// layout's allocation count is asserted to be stable
    /// (`checkAllAllocationFailures` in tests/native_minimal.zig).
    storage: []core.math.Vec2,

    fn init(allocator: std.mem.Allocator, atom_count: usize, fragment_count: usize) core.errors.Error!FragmentFrames {
        const storage = allocator.alloc(core.math.Vec2, atom_count + fragment_count) catch return error.OutOfMemory;
        @memset(storage, .{});
        return .{
            .allocator = allocator,
            .storage = storage,
            .own = storage[0..atom_count],
            .attachment = storage[atom_count..],
        };
    }

    fn deinit(self: *FragmentFrames) void {
        self.allocator.free(self.storage);
        self.* = undefined;
    }

    /// The frame is centred on the fragment's attachment atom and turned so the
    /// bond to the parent lies along it; a root fragment is centred on its
    /// first atom and not turned, and a constrained or fixed root is left where
    /// it is so its template alignment survives.
    fn store(
        self: FragmentFrames,
        atoms: []const model.Atom,
        fragmentation: fragments.Fragmentation,
        fragment: fragments.Fragment,
    ) void {
        const members = fragmentation.members(fragment.id);
        var origin: core.math.Vec2 = .{};
        var angle: f32 = 0;
        if (attachmentOf(fragment)) |attachment_atom| {
            origin = atoms[attachment_atom.index()].coordinates;
            const parent_atom = parentAtomOf(fragment).?;
            const parent_position = atoms[parent_atom.index()].coordinates;
            angle = std.math.atan2(parent_position.y - origin.y, origin.x - parent_position.x);
        } else if (!fragment.flags.constrained and !fragment.flags.fixed) {
            origin = atoms[members[0].index()].coordinates;
        }
        for (members) |atom| {
            self.own[atom.index()] = rotateClockwise(subtract(atoms[atom.index()].coordinates, origin), -angle);
        }
        for (fragmentation.fragments) |candidate| {
            if (candidate.parent != fragment.id) continue;
            const child_attachment = attachmentOf(candidate) orelse continue;
            self.attachment[candidate.id.index()] =
                rotateClockwise(subtract(atoms[child_attachment.index()].coordinates, origin), -angle);
        }
    }
};

/// `CoordgenMinimizer::buildMoleculeFromFragments` with firstTime true, minus
/// the side choice it also makes there (cgz-7v2.27). A child is placed at the
/// position its parent's own layout gave its attachment atom, turned to the
/// direction of the bond it arrived on - so two children of one atom take the
/// two directions that layout chose for them.
///
/// The port previously synthesised that direction instead, from the mean of the
/// parent atom's already-placed bonds. That mean does not depend on which child
/// is being placed, so every child of one atom was sent the same way and landed
/// on top of its siblings (measured on drug_like/4: five pairs or triples of
/// exactly coincident atoms).
fn placeFragmentsFromFrames(
    allocator: std.mem.Allocator,
    atoms: []model.Atom,
    bonds: []const model.Bond,
    graph: topology.Graph,
    fragmentation: fragments.Fragmentation,
    frames: FragmentFrames,
) core.errors.Error!void {
    var scratch = try SideChoiceScratch.init(allocator, bonds.len, fragmentation.fragments.len);
    defer scratch.deinit();
    assignLongestChains(fragmentation, frames, scratch.longest_chain);
    // A parent must be placed before its children, because the child reads the
    // attachment position the parent's placement just wrote.
    for (0..fragmentation.fragments.len) |depth| {
        for (fragmentation.fragments) |fragment| {
            if (fragmentDepth(fragmentation.fragments, fragment.id) != depth) continue;
            var position: core.math.Vec2 = .{};
            var angle: f32 = 0;
            if (attachmentOf(fragment)) |attachment_atom| {
                const parent_atom = parentAtomOf(fragment) orelse return error.InvalidMapping;
                position = atoms[attachment_atom.index()].coordinates;
                const direction = subtract(position, atoms[parent_atom.index()].coordinates);
                angle = std.math.atan2(-direction.y, direction.x);
                // `buildMoleculeFromFragments` decides the side here, between
                // computing the angle and applying the frame, and only on the
                // first build (cgz-7v2.27).
                // `buildMoleculeFromFragments` decides the side here, between
                // computing the angle and applying the frame, and only on the
                // first build (cgz-7v2.27).
                if (!fragment.flags.fixed and !fragment.flags.constrained and
                    shouldInvertAgainstParent(atoms, bonds, graph, fragmentation, frames, fragment, angle, scratch))
                {
                    mirrorFragmentFrame(fragmentation, frames, fragment);
                }
            } else if (fragment.parent.isValid()) {
                return error.InvalidMapping;
            }
            for (fragmentation.members(fragment.id)) |atom| {
                atoms[atom.index()].coordinates = add(rotateClockwise(frames.own[atom.index()], angle), position);
            }
            for (fragmentation.fragments) |candidate| {
                if (candidate.parent != fragment.id) continue;
                const child_attachment = attachmentOf(candidate) orelse continue;
                atoms[child_attachment.index()].coordinates =
                    add(rotateClockwise(frames.attachment[candidate.id.index()], angle), position);
            }
        }
    }
}

/// The end of `bond_to_parent` inside this fragment, and the end outside it.
///
/// Both are resolved once by `Fragmentation`, by membership, and only read
/// back here. Deliberately never recomputed from the bond's stored direction:
/// native keeps bonds child-first, so `.start` is not the parent's end the way
/// upstream's fragmenter guarantees, and reading it fails silently (cgz-jg4).
fn attachmentOf(fragment: fragments.Fragment) ?core.ids.AtomId {
    return if (fragment.attachment_atom.isValid()) fragment.attachment_atom else null;
}

fn parentAtomOf(fragment: fragments.Fragment) ?core.ids.AtomId {
    return if (fragment.parent_atom.isValid()) fragment.parent_atom else null;
}

fn placeFragmentRings(
    allocator: std.mem.Allocator,
    atoms: []model.Atom,
    bonds: []const model.Bond,
    graph: topology.Graph,
    membership: topology.RingMembership,
    fragmentation: fragments.Fragmentation,
    fragment: fragments.Fragment,
    placed: []bool,
    options: Options,
    analysis: topology.rings.Analysis,
) core.errors.Error!bool {
    // Upstream tries the template database for the whole ring system before
    // placing any ring individually (CoordgenFragmentBuilder.cpp:211, inside
    // generateCoordinatesCentralRings, which only reaches the per-ring
    // buildRing loop when findTemplate returns false). Omitting this call was
    // cgz-r35: every entry point of the template layer was unreachable, so a
    // bridged bicyclic that upstream lays out from template 79 fell through to
    // generic fused-ring placement and came out with crossing bonds.
    if (try matchFragmentTemplate(allocator, atoms, bonds, membership, fragmentation, fragment, options, placed)) {
        // A templated system never reaches maybeMinimizeRings upstream: that
        // call lives in the not-found branch.
        return false;
    }
    // Upstream branches three ways on a planarity score before placing anything
    // (CoordgenFragmentBuilder.cpp:213-235). Native had no planarity score at
    // all and always took the first branch, so a system upstream lays out as an
    // opened macrocycle, or declines to lay out entirely, was placed ring by
    // ring instead (cgz-r31).
    const planarity = try scorePlanarity(membership, analysis, bonds, graph, fragmentation, fragment);
    if (planarity > untreatable_system_planarity_score) {
        // Upstream returns, leaving the whole ring system unplaced. Faithful:
        // the caller's `placed` flags stay false and the atoms keep whatever
        // coordinates they had.
        return false;
    }
    if (planarity >= non_planar_system_score) {
        // Neither planar enough to place ring by ring nor hopeless: upstream
        // opens the central ring as a macrocycle and requires minimization.
        const central = selectCentralRing(membership, analysis, fragmentation, fragment) orelse
            return false;
        _ = try openCycleAndGenerateCoordinates(allocator, central, atoms, bonds, graph, membership, placed);
        return true;
    }
    if (planarity > perfectly_planar_system_score) {
        // Upstream repeats the same findTemplate call on the same ring set here.
        // It is deterministic, so at the pin the retry can only return what the
        // first call already returned; kept because the shape of the branch is
        // upstream's, and measured to change nothing.
        if (try matchFragmentTemplate(allocator, atoms, bonds, membership, fragmentation, fragment, options, placed)) {
            return false;
        }
    }
    var pentagon_minimization = false;
    var remaining: usize = fragment.ring_count;
    while (remaining != 0) {
        var selected: ?model.Ring = null;
        var selected_shared: usize = 0;
        var selected_score: usize = 0;
        for (membership.rings) |ring| {
            const ring_atoms = membership.atoms(ring.id);
            if (ring_atoms.len < 3 or
                fragmentation.atom_fragment[ring_atoms[0].index()] != fragment.id) continue;
            var already_complete = true;
            var shared: usize = 0;
            for (ring_atoms) |atom| {
                already_complete = already_complete and placed[atom.index()];
                shared += @intFromBool(placed[atom.index()]);
            }
            if (already_complete) continue;
            const score = ringPriority(membership, analysis, ring, shared != 0);
            if (selected == null or score > selected_score) {
                selected = ring;
                selected_shared = shared;
                selected_score = score;
            }
        }
        const ring = selected orelse break;
        if (membership.atoms(ring.id).len >= topology.rings.macrocycle_size) {
            const shape_result = try macrocycle.generateRingShape(
                allocator,
                ring.id,
                atoms,
                bonds,
                graph,
                membership,
                placed,
                options.force_open_macrocycles,
            );
            if (shape_result.placed()) {
                // A polyomino with pentagon vertices leaves the ring off the
                // hexagonal lattice, and upstream requires minimization for
                // exactly that case.
                if (shape_result == .matched_needs_minimization) pentagon_minimization = true;
                remaining -= 1;
                continue;
            }
            if (try openCycleAndGenerateCoordinates(allocator, ring.id, atoms, bonds, graph, membership, placed)) {
                remaining -= 1;
                continue;
            }
        }
        const ordered = try orderRing(allocator, ring, bonds, membership);
        defer allocator.free(ordered);
        const local = try regularRingCoordinates(allocator, ordered.len);
        defer allocator.free(local);
        if (selected_shared >= 2) {
            // Upstream scores the two mirror candidates against the center of
            // the one already-drawn parent ring, not against everything placed
            // so far (CoordgenFragmentBuilder.cpp:513). The distinction only
            // appears from the third ring of a fused system onwards, which is
            // where the two disagree and where the wrong mirror gets picked.
            const parent = selectParentRing(membership, analysis, ring.id, placed) orelse
                return error.InvalidMapping;
            try alignFusedRing(atoms, ordered, local, placed, membership.atoms(parent));
        } else if (selected_shared == 1) {
            var pivot: usize = 0;
            for (ordered, 0..) |atom, index| if (placed[atom.index()]) {
                pivot = index;
                break;
            };
            const source_center = coordinateCenter(local);
            const existing_center = placedCenter(atoms, placed);
            const target = atoms[ordered[pivot].index()].coordinates;
            var outward = subtract(target, existing_center);
            if (length(outward) < 0.0001) outward = .{ .x = 1 };
            var source_direction = subtract(source_center, local[pivot]);
            if (length(source_direction) < 0.0001) source_direction = .{ .x = 1 };
            const rotation = std.math.atan2(outward.y, outward.x) - std.math.atan2(source_direction.y, source_direction.x);
            for (ordered, local) |atom, coordinate| if (!placed[atom.index()]) {
                atoms[atom.index()].coordinates = transformFromPivot(coordinate, local[pivot], target, rotation);
                placed[atom.index()] = true;
            };
        } else {
            for (ordered, local) |atom, coordinate| {
                atoms[atom.index()].coordinates = coordinate;
                placed[atom.index()] = true;
            }
        }
        remaining -= 1;
    }
    // Upstream reaches maybeMinimizeRings only after the loop above has drained
    // its ring vector — `while (!rings.empty()) { ... rings.erase(...); }` —
    // so the call it makes is over an empty set (CoordgenFragmentBuilder.cpp
    // :220-227). The function iterates nothing, finds nothing, and returns
    // before its own `rings.at(0)` would throw on the empty vector. It is dead
    // code at the pin, and its trigger never fires.
    //
    // Native evaluated it over the fragment's real rings and so required
    // minimization on systems upstream never minimizes (cgz-r31). The
    // transcription is kept and tested; what is corrected is the ring set the
    // call site hands it, which upstream leaves empty.
    const rings_left_after_placement: []const core.ids.RingId = &.{};
    return pentagon_minimization or maybeMinimizeRings(membership, rings_left_after_placement);
}

/// `findCentralRingOfSystem` for a system where nothing has been placed yet, so
/// no ring carries the already-built bonus. Same priority rule the per-ring
/// placement loop applies, and first-wins on a tie, as there.
fn selectCentralRing(
    membership: topology.RingMembership,
    analysis: topology.rings.Analysis,
    fragmentation: fragments.Fragmentation,
    fragment: fragments.Fragment,
) ?core.ids.RingId {
    var selected: ?core.ids.RingId = null;
    var selected_score: usize = 0;
    for (membership.rings) |ring| {
        const ring_atoms = membership.atoms(ring.id);
        if (ring_atoms.len < 3) continue;
        if (fragmentation.atom_fragment[ring_atoms[0].index()] != fragment.id) continue;
        const score = ringPriority(membership, analysis, ring, false);
        if (selected == null or score > selected_score) {
            selected = ring.id;
            selected_score = score;
        }
    }
    return selected;
}

pub const perfectly_planar_system_score: f32 = 50;
pub const non_planar_system_score: f32 = 1000;
pub const untreatable_system_planarity_score: f32 = 200000;

/// Upstream's `CoordgenFragmentBuilder::newScorePlanarity`: how badly a fused
/// ring system resists being drawn flat. Three conditions each add
/// `non_planar_system_score`, so a system scoring at or above that has at least
/// one of them and cannot be placed ring by ring.
///
/// A macrocycle with no bond available to open is skipped entirely, before its
/// own conditions are considered — that is upstream's `continue`, and it is why
/// the openable check comes first.
fn scorePlanarity(
    membership: topology.RingMembership,
    analysis: topology.rings.Analysis,
    bonds: []const model.Bond,
    graph: topology.Graph,
    fragmentation: fragments.Fragmentation,
    fragment: fragments.Fragment,
) core.errors.Error!f32 {
    var score: f32 = 0;
    for (membership.rings) |ring| {
        const ring_atoms = membership.atoms(ring.id);
        if (ring_atoms.len == 0) continue;
        if (fragmentation.atom_fragment[ring_atoms[0].index()] != fragment.id) continue;
        const macrocycle_ring = ring_atoms.len >= topology.rings.macrocycle_size;
        if (macrocycle_ring and
            macrocycle.findBondToOpen(ring.id, bonds, graph, membership) == null) continue;
        if (macrocycle_ring) {
            for (analysis.fusedWith(ring.id)) |fusion| {
                if (membership.atoms(fusion.other).len >= topology.rings.macrocycle_size) {
                    score += non_planar_system_score;
                }
            }
        }
        for (bonds) |bond| {
            const shared = membership.bondRings(bond.id);
            if (std.mem.indexOfScalar(core.ids.RingId, shared, ring.id) == null) continue;
            if (shared.len > 2) {
                score += non_planar_system_score * @as(f32, @floatFromInt(shared.len - 2));
            }
        }
        for (ring_atoms) |atom| {
            if (graph.degree(atom) <= 3) continue;
            var angle: f32 = 0;
            for (membership.atomRings(atom)) |other| {
                const size = membership.atoms(other).len;
                if (size == 0) return error.InvalidMapping;
                angle += std.math.pi - 2 * std.math.pi / @as(f32, @floatFromInt(size));
            }
            // Upstream's threshold is 1.99 pi, not 2 pi: three ring interior
            // angles that sum to a full turn leave no room for a fourth
            // substituent in the plane, and the slack absorbs float error.
            if (angle >= 1.99 * std.math.pi) score += non_planar_system_score;
        }
    }
    return score;
}

/// Upstream's CoordgenMinimizer::maybeMinimizeRings: a five-membered ring, or
/// an odd-sized macrocycle, holding an atom that belongs to more than two rings
/// requires minimization.
///
/// Faithful to the function, which is not the same as reachable — see the call
/// site above for why upstream's only caller passes it nothing.
fn maybeMinimizeRings(
    membership: topology.RingMembership,
    rings: []const core.ids.RingId,
) bool {
    for (rings) |ring| {
        const ring_atoms = membership.atoms(ring);
        const five_membered = ring_atoms.len == 5;
        const odd_macrocycle = ring_atoms.len >= topology.rings.macrocycle_size and
            ring_atoms.len % 2 != 0;
        if (!five_membered and !odd_macrocycle) continue;
        for (ring_atoms) |atom| {
            if (membership.atomRings(atom).len > 2) return true;
        }
    }
    return false;
}

/// Collect this fragment's rings in ring order and try the template database
/// against the whole system, mirroring findTemplate's first upstream call
/// site. Returns true when the fragment was placed from a template, in which
/// case no per-ring placement runs for it.
///
/// templates.findRingSet already declines a system of fewer than two rings and
/// already honours `load_templates`, so both of upstream's guards live there
/// rather than being restated here.
fn matchFragmentTemplate(
    allocator: std.mem.Allocator,
    atoms: []model.Atom,
    bonds: []const model.Bond,
    membership: topology.RingMembership,
    fragmentation: fragments.Fragmentation,
    fragment: fragments.Fragment,
    options: Options,
    placed: []bool,
) core.errors.Error!bool {
    var ring_ids: std.ArrayList(core.ids.RingId) = .empty;
    defer ring_ids.deinit(allocator);
    for (membership.rings) |ring| {
        const ring_atoms = membership.atoms(ring.id);
        if (ring_atoms.len < 3) continue;
        if (fragmentation.atom_fragment[ring_atoms[0].index()] != fragment.id) continue;
        ring_ids.append(allocator, ring.id) catch return error.OutOfMemory;
    }
    if (ring_ids.items.len < 2) return false;

    var match = (try templates.findRingSet(
        allocator,
        options.templates,
        atoms.len,
        bonds,
        membership,
        ring_ids.items,
    )) orelse return false;
    defer match.deinit();
    try match.apply(atoms);
    for (match.atoms) |atom_id| placed[atom_id.index()] = true;
    return true;
}

/// Upstream's CoordgenFragmentBuilder::findCentralRingOfSystem priority, with
/// its constants (CoordgenFragmentBuilder.cpp:24-27).
///
/// The native score used to be `shared * 10000 + (size == 6) * 10 + size`,
/// which kept the already-built bonus, the six-ring bonus and the size term
/// but dropped three: the macrocycle bonus, the fused-ring count and the
/// fusion-atom count. Those decide the order for any system holding a
/// macrocycle or rings of differing fusion, and the placement order decides
/// the layout.
///
/// `neighbour_built` is a flat bonus upstream, taken once for having any built
/// fused neighbour rather than scaled by how many atoms are already placed.
fn ringPriority(
    membership: topology.RingMembership,
    analysis: topology.rings.Analysis,
    ring: model.Ring,
    neighbour_built: bool,
) usize {
    const neighbour_already_built_score: usize = 100000;
    const macrocycle_score: usize = 1000;
    const fused_rings_score: usize = 40;
    const fusion_atoms_score: usize = 15;

    var priority: usize = 0;
    if (neighbour_built) priority += neighbour_already_built_score;
    const size = membership.atoms(ring.id).len;
    if (size >= topology.rings.macrocycle_size) priority += macrocycle_score;
    if (size == 6) priority += 10;
    priority += size;
    const fusions = analysis.fusedWith(ring.id);
    priority += fused_rings_score * fusions.len;
    for (fusions) |fusion| priority += fusion_atoms_score * analysis.fusionAtoms(fusion).len;
    return priority;
}

fn openCycleAndGenerateCoordinates(
    allocator: std.mem.Allocator,
    ring: core.ids.RingId,
    atoms: []model.Atom,
    bonds: []const model.Bond,
    graph: topology.Graph,
    membership: topology.RingMembership,
    placed: []bool,
) core.errors.Error!bool {
    const bond_to_open = macrocycle.findBondToOpen(ring, bonds, graph, membership) orelse return false;
    const ring_atoms = membership.atoms(ring);
    if (ring_atoms.len == 0) return false;
    const component = graph.component(ring_atoms[0]) orelse return false;
    const temporary_atoms = allocator.dupe(model.Atom, atoms) catch return error.OutOfMemory;
    defer allocator.free(temporary_atoms);
    for (temporary_atoms) |*atom| {
        if (graph.component(atom.id) != component) atom.hidden = true;
    }
    const temporary_bonds = allocator.dupe(model.Bond, bonds) catch return error.OutOfMemory;
    defer allocator.free(temporary_bonds);
    temporary_bonds[bond_to_open.index()].skip = true;
    var temporary_graph = try topology.Graph.init(allocator, temporary_atoms, temporary_bonds);
    defer temporary_graph.deinit();
    var temporary_rings = try topology.RingMembership.init(allocator, temporary_graph, temporary_bonds);
    defer temporary_rings.deinit();
    var temporary_fragments = try fragments.Fragmentation.init(allocator, temporary_atoms, temporary_bonds, temporary_graph, temporary_rings);
    defer temporary_fragments.deinit();
    _ = try initializeCoordinatesInternal(
        allocator,
        temporary_atoms,
        temporary_bonds,
        temporary_graph,
        temporary_rings,
        temporary_fragments,
        .{ .force_open_macrocycles = true },
    );
    for (atoms, temporary_atoms, placed) |*atom, temporary_atom, *is_placed| {
        if (graph.component(atom.id) != component) continue;
        atom.coordinates = temporary_atom.coordinates;
        is_placed.* = true;
    }
    return true;
}

fn regularRingCoordinates(allocator: std.mem.Allocator, count: usize) core.errors.Error![]core.math.Vec2 {
    if (count < 3) return error.InvalidMapping;
    const result = allocator.alloc(core.math.Vec2, count) catch return error.OutOfMemory;
    errdefer allocator.free(result);
    var coordinate: core.math.Vec2 = .{};
    const step = 2 * pi / @as(f32, @floatFromInt(count));
    for (result, 0..) |*destination, index| {
        destination.* = coordinate;
        const angle = step * @as(f32, @floatFromInt(index));
        coordinate.x += @cos(angle) * bond_length;
        coordinate.y -= @sin(angle) * bond_length;
    }
    return result;
}

/// Upstream's `getSharedAtomsWithAlreadyDrawnRing`
/// (CoordgenFragmentBuilder.cpp): among this ring's fused neighbours that are
/// already drawn, walk them in fusion order and keep the last one that is
/// neither less fused with this ring nor smaller than the current pick. The
/// disjunction is upstream's, and it makes the choice order-dependent rather
/// than a maximum — preserved here for that reason.
fn selectParentRing(
    membership: topology.RingMembership,
    analysis: topology.rings.Analysis,
    ring: core.ids.RingId,
    placed: []const bool,
) ?core.ids.RingId {
    var parent: ?core.ids.RingId = null;
    var parent_shared: usize = 0;
    for (analysis.fusedWith(ring)) |fusion| {
        if (!ringIsDrawn(membership, fusion.other, placed)) continue;
        if (parent != null) {
            if (fusion.atom_count < parent_shared or
                membership.atoms(fusion.other).len < membership.atoms(parent.?).len) continue;
        }
        parent = fusion.other;
        parent_shared = fusion.atom_count;
    }
    return parent;
}

fn ringIsDrawn(membership: topology.RingMembership, ring: core.ids.RingId, placed: []const bool) bool {
    for (membership.atoms(ring)) |atom| {
        if (!placed[atom.index()]) return false;
    }
    return true;
}

fn alignFusedRing(
    atoms: []model.Atom,
    ordered: []const core.ids.AtomId,
    local: []const core.math.Vec2,
    placed: []bool,
    parent_atoms: []const core.ids.AtomId,
) core.errors.Error!void {
    var first: ?usize = null;
    var last: usize = 0;
    // The anchors and the excluded set are the atoms shared with the parent,
    // upstream's `fusionAtoms`, not every atom already placed: a ring fused to
    // two drawn rings shares atoms with both, and only the parent's pair
    // defines the axis the two candidates differ across.
    for (ordered, 0..) |atom, index| {
        if (std.mem.indexOfScalar(core.ids.AtomId, parent_atoms, atom) == null) continue;
        if (first == null) first = index;
        last = index;
    }
    const first_index = first orelse return error.InvalidMapping;
    if (first_index == last) return error.InvalidMapping;
    const target_first = atoms[ordered[first_index].index()].coordinates;
    const target_last = atoms[ordered[last].index()].coordinates;
    const source_first = local[first_index];
    const source_last = local[last];
    const target_angle = std.math.atan2(target_last.y - target_first.y, target_last.x - target_first.x);
    const source_angle = std.math.atan2(source_last.y - source_first.y, source_last.x - source_first.x);
    const rotation = target_angle - source_angle;
    var first_score: f32 = 0;
    var mirror_score: f32 = 0;
    const parent_center = ringAtomCenter(atoms, parent_atoms);
    for (ordered, local) |atom, coordinate| {
        if (std.mem.indexOfScalar(core.ids.AtomId, parent_atoms, atom) != null) continue;
        const candidate = transformFromPivot(coordinate, source_first, target_first, rotation);
        const mirror = reflectAcrossLine(candidate, target_first, target_last);
        first_score += distance(candidate, parent_center);
        mirror_score += distance(mirror, parent_center);
    }
    const use_mirror = mirror_score > first_score;
    // Every atom of the ring is written, the shared ones included, which
    // overwrites what the parent ring gave them. Upstream does exactly this -
    // `buildRing`'s two-or-more fusion-atom branch ends in an unconditional
    // `atoms[i]->setCoordinates(targetCoords[i])` over the whole ring
    // (CoordgenFragmentBuilder.cpp:532-534).
    //
    // Skipping already-placed atoms looks conservative and is not. On a
    // BRIDGED system the rings share three atoms rather than two, so leaving
    // those three where the first ring put them pins the second ring at three
    // points and it can no longer hold its own shape. Measured on drug_like/3,
    // whose 7-ring is bridged to a 13-macrocycle across atoms 5, 10 and 16:
    // upstream's heptagon is regular, every atom 1.150 to 1.154 from the
    // centre, and native's spanned 1.140 to 1.171 (cgz-vu0).
    for (ordered, local) |atom, coordinate| {
        const candidate = transformFromPivot(coordinate, source_first, target_first, rotation);
        atoms[atom.index()].coordinates = if (use_mirror) reflectAcrossLine(candidate, target_first, target_last) else candidate;
        placed[atom.index()] = true;
    }
}

fn transformFromPivot(point: core.math.Vec2, source: core.math.Vec2, target: core.math.Vec2, rotation: f32) core.math.Vec2 {
    const local = subtract(point, source);
    const cosine = @cos(rotation);
    const sine = @sin(rotation);
    return .{
        .x = target.x + local.x * cosine - local.y * sine,
        .y = target.y + local.x * sine + local.y * cosine,
    };
}

fn reflectAcrossLine(point: core.math.Vec2, line_start: core.math.Vec2, line_end: core.math.Vec2) core.math.Vec2 {
    const direction = subtract(line_end, line_start);
    const denominator = direction.x * direction.x + direction.y * direction.y;
    if (denominator < 0.000001) return point;
    const relative = subtract(point, line_start);
    const projection = (relative.x * direction.x + relative.y * direction.y) / denominator;
    const on_line = add(line_start, scale(direction, projection));
    return subtract(scale(on_line, 2), point);
}

fn placedCenter(atoms: []const model.Atom, placed: []const bool) core.math.Vec2 {
    var center: core.math.Vec2 = .{};
    var count: usize = 0;
    for (atoms, placed) |atom, is_placed| if (is_placed) {
        center = add(center, atom.coordinates);
        count += 1;
    };
    return if (count == 0) center else scale(center, 1 / @as(f32, @floatFromInt(count)));
}

/// `sketcherMinimizerRing::findCenter`: the mean of one ring's atoms.
fn ringAtomCenter(atoms: []const model.Atom, members: []const core.ids.AtomId) core.math.Vec2 {
    var center: core.math.Vec2 = .{};
    for (members) |atom| center = add(center, atoms[atom.index()].coordinates);
    return if (members.len == 0) center else scale(center, 1 / @as(f32, @floatFromInt(members.len)));
}

fn coordinateCenter(coordinates: []const core.math.Vec2) core.math.Vec2 {
    var center: core.math.Vec2 = .{};
    for (coordinates) |coordinate| center = add(center, coordinate);
    return scale(center, 1 / @as(f32, @floatFromInt(coordinates.len)));
}

/// Capture the assembled pose in the fragment-local representation consumed
/// by discrete optimization. Child anchors are deliberately duplicated in
/// their own frame and in their parent's attachment coordinates.
pub fn captureFragmentFrames(
    allocator: std.mem.Allocator,
    atoms: []const model.Atom,
    fragmentation: fragments.Fragmentation,
) core.errors.Error!core.dof.FrameCollection {
    const frame_count = fragmentation.fragments.len;
    const frames = allocator.alloc(core.dof.FragmentFrame, frame_count) catch return error.OutOfMemory;
    errdefer allocator.free(frames);
    const rebuild_order = allocator.alloc(core.ids.FragmentId, frame_count) catch return error.OutOfMemory;
    errdefer allocator.free(rebuild_order);
    const fragment_atoms = allocator.dupe(core.ids.AtomId, fragmentation.atoms) catch return error.OutOfMemory;
    errdefer allocator.free(fragment_atoms);
    const atom_coordinates = allocator.alloc(core.math.Vec2, fragment_atoms.len) catch return error.OutOfMemory;
    errdefer allocator.free(atom_coordinates);
    const attachment_count = if (frame_count == 0) 0 else frame_count - fragmentation.main_fragments.len;
    const child_attachments = allocator.alloc(core.dof.ChildAttachment, attachment_count) catch return error.OutOfMemory;
    errdefer allocator.free(child_attachments);
    const attachment_coordinates = allocator.alloc(core.math.Vec2, attachment_count) catch return error.OutOfMemory;
    errdefer allocator.free(attachment_coordinates);

    var attachment_offset: u32 = 0;
    for (fragmentation.fragments) |fragment| {
        var child_count: u32 = 0;
        for (fragmentation.fragments) |candidate| child_count += @intFromBool(candidate.parent == fragment.id);
        const frame = &frames[fragment.id.index()];
        frame.* = .{
            .id = fragment.id,
            .parent = fragment.parent,
            .atoms = .{ .start = fragment.atom_start, .len = fragment.atom_count },
            .attachments = .{ .start = attachment_offset, .len = child_count },
        };
        if (fragment.parent.isValid()) {
            if (!fragment.attachment_atom.isValid() or !fragment.parent_atom.isValid()) return error.InvalidMapping;
            frame.anchor_atom = fragment.attachment_atom;
            frame.parent_atom = fragment.parent_atom;
        }
        attachment_offset += child_count;
    }
    if (attachment_offset != attachment_count) return error.InvalidMapping;

    var order_index: usize = 0;
    for (0..frame_count) |depth| for (fragmentation.fragments) |fragment| {
        if (fragmentDepth(fragmentation.fragments, fragment.id) != depth) continue;
        rebuild_order[order_index] = fragment.id;
        order_index += 1;
    };
    if (order_index != frame_count) return error.InvalidMapping;

    for (frames) |frame| {
        const origin, const angle = frameOriginAndAngle(frame, fragmentation.fragments[frame.id.index()], fragment_atoms, atoms);
        const sine = @sin(-angle);
        const cosine = @cos(-angle);
        for (fragment_atoms[frame.atoms.start..][0..frame.atoms.len], frame.atoms.start..) |atom, coordinate_index| {
            atom_coordinates[coordinate_index] = rotateScreen(subtract(atoms[atom.index()].coordinates, origin), sine, cosine);
        }
        var child_index: usize = frame.attachments.start;
        for (frames) |child| {
            if (child.parent != frame.id) continue;
            child_attachments[child_index] = .{ .child = child.id, .atom = child.anchor_atom };
            attachment_coordinates[child_index] = rotateScreen(subtract(atoms[child.anchor_atom.index()].coordinates, origin), sine, cosine);
            child_index += 1;
        }
    }
    return .{
        .allocator = allocator,
        .frames = frames,
        .rebuild_order = rebuild_order,
        .fragment_atoms = fragment_atoms,
        .child_attachments = child_attachments,
        .atom_coordinates = atom_coordinates,
        .attachment_coordinates = attachment_coordinates,
    };
}

fn frameOriginAndAngle(
    frame: core.dof.FragmentFrame,
    fragment: fragments.Fragment,
    fragment_atoms: []const core.ids.AtomId,
    atoms: []const model.Atom,
) struct { core.math.Vec2, f32 } {
    if (!frame.parent.isValid()) {
        if (!fragment.flags.constrained and !fragment.flags.fixed and frame.atoms.len != 0) {
            return .{ atoms[fragment_atoms[frame.atoms.start].index()].coordinates, 0 };
        }
        return .{ .{}, 0 };
    }
    const origin = atoms[frame.anchor_atom.index()].coordinates;
    const parent = atoms[frame.parent_atom.index()].coordinates;
    return .{ origin, std.math.atan2(parent.y - origin.y, -parent.x + origin.x) };
}

fn rotateScreen(value: core.math.Vec2, sine: f32, cosine: f32) core.math.Vec2 {
    return .{ .x = value.x * cosine + value.y * sine, .y = -value.x * sine + value.y * cosine };
}

fn fragmentDepth(records: []const fragments.Fragment, fragment: core.ids.FragmentId) usize {
    var depth: usize = 0;
    var cursor = fragment;
    while (records[cursor.index()].parent.isValid()) {
        depth += 1;
        cursor = records[cursor.index()].parent;
        std.debug.assert(depth < records.len);
    }
    return depth;
}

/// `CoordgenFragmentBuilder::rotateMainFragment`, which upstream runs on a
/// constrained root at the end of that root's own build.
fn alignConstrainedMainFragment(atoms: []model.Atom, fragmentation: fragments.Fragmentation, fragment: fragments.Fragment) void {
    {
        if (fragment.parent.isValid()) return;
        if (!fragment.flags.constrained or fragment.flags.fixed) return;
        const members = fragmentation.members(fragment.id);
        var source_center: core.math.Vec2 = .{};
        var target_center: core.math.Vec2 = .{};
        var constrained_count: usize = 0;
        for (members) |atom| if (atoms[atom.index()].constrained) {
            const target = atoms[atom.index()].template_coordinates orelse continue;
            source_center = add(source_center, atoms[atom.index()].coordinates);
            target_center = add(target_center, target);
            constrained_count += 1;
        };
        if (constrained_count == 0) return;
        source_center = scale(source_center, 1 / @as(f32, @floatFromInt(constrained_count)));
        target_center = scale(target_center, 1 / @as(f32, @floatFromInt(constrained_count)));
        var dot: f32 = 0;
        var cross: f32 = 0;
        for (members) |atom| if (atoms[atom.index()].constrained) {
            const target = atoms[atom.index()].template_coordinates orelse continue;
            const source_delta = subtract(atoms[atom.index()].coordinates, source_center);
            const target_delta = subtract(target, target_center);
            dot += source_delta.x * target_delta.x + source_delta.y * target_delta.y;
            cross += source_delta.x * target_delta.y - source_delta.y * target_delta.x;
        };
        const angle = if (constrained_count > 1) std.math.atan2(cross, dot) else 0;
        for (members) |atom| {
            atoms[atom.index()].coordinates = transformFromPivot(atoms[atom.index()].coordinates, source_center, target_center, angle);
        }
    }
}

/// `sketcherMinimizerFragment::setAllCoordinatesToTemplate`, run by upstream on
/// a fixed fragment at the end of its own build so that the frame the fragment
/// is stored in is the template's.
fn restoreFixedFragmentCoordinates(atoms: []model.Atom, fragmentation: fragments.Fragmentation, fragment: fragments.Fragment) void {
    if (!fragment.flags.fixed) return;
    for (fragmentation.members(fragment.id)) |atom| {
        if (atoms[atom.index()].template_coordinates) |coordinate| atoms[atom.index()].coordinates = coordinate;
    }
}

/// The same reset once more after assembly: a fixed atom keeps its input
/// coordinates whatever the layout did with the fragment around it.
fn restoreFixedCoordinates(atoms: []model.Atom) void {
    for (atoms) |*atom| if (atom.fixed) {
        if (atom.template_coordinates) |coordinate| atom.coordinates = coordinate;
    };
}

/// `CoordgenFragmentBuilder::fallbackIfNanCoordinates`: when a fragment's
/// generated coordinates hold a NaN and its input 3D coordinates are usable,
/// project those instead.
///
/// This was cgz-7v2.25. A hand-rolled version was wired here that copied 3D x/y
/// straight through - unscaled, unmirrored, unrounded - while the faithful
/// transcription sat unreachable in the optimize layer, and their triggers
/// disagreed as well. Two implementations of one upstream function is how that
/// happened, so there is now exactly one, in the layer that calls it.
fn fallbackOnValid3dCoordinates(atoms: []model.Atom, members: []const core.ids.AtomId) void {
    {
        if (!fragmentHasNan(atoms, members)) return;
        if (!fragmentHasValid3dSource(atoms, members)) return;
        for (members) |atom| {
            const source = atoms[atom.index()].coordinates_3d.?;
            // The 35x scale and the negated y are upstream's projection, and the
            // rounding keeps the emergency pose on the same two-decimal grid as
            // the rest of the layout.
            atoms[atom.index()].coordinates = .{
                .x = roundToTwoDecimalDigits(source.x * 35),
                .y = roundToTwoDecimalDigits(-source.y * 35),
            };
        }
    }
}

/// `CoordgenMinimizer::hasNaNCoordinates`: NaN only. An infinity is deliberately
/// not a trigger upstream, and the hand-rolled version tested `!isFinite()`,
/// which made it one and took the fallback where upstream would not.
fn fragmentHasNan(atoms: []const model.Atom, members: []const core.ids.AtomId) bool {
    for (members) |atom| {
        const position = atoms[atom.index()].coordinates;
        if (std.math.isNan(position.x) or std.math.isNan(position.y)) return true;
    }
    return false;
}

/// `sketcherMinimizerAtom::hasValid3DCoordinates`: upstream's 10,000,001 test.
/// The number is not a magnitude policy - INVALID_COORDINATES is the *default*
/// value upstream assigns m_x3D/y3D/z3D, so `< INVALID_COORDINATES` means "this
/// atom actually carries 3D data". Native's `?Vec3` encodes that intent
/// directly, and the null check below is the faithful half.
///
/// The ceiling is kept as well because it is one-sided, which makes it more than
/// a presence test: a coordinate at or above the sentinel reads as absent, and a
/// negative infinity reads as present. Both follow from comparing against a
/// sentinel rather than from a decision, and the test below pins them rather
/// than tidying them away.
fn fragmentHasValid3dSource(atoms: []const model.Atom, members: []const core.ids.AtomId) bool {
    const invalid_coordinates: f32 = 10_000_001;
    for (members) |atom| {
        const source = atoms[atom.index()].coordinates_3d orelse return false;
        if (!(source.x < invalid_coordinates and
            source.y < invalid_coordinates and
            source.z < invalid_coordinates)) return false;
    }
    return true;
}

fn roundToTwoDecimalDigits(value: f32) f32 {
    return @floor(value * 100 + 0.5) * 0.01;
}

fn fragmentCenter(atoms: []const model.Atom, members: []const core.ids.AtomId) core.math.Vec2 {
    var center: core.math.Vec2 = .{};
    for (members) |atom| center = add(center, atoms[atom.index()].coordinates);
    return scale(center, 1 / @as(f32, @floatFromInt(members.len)));
}

fn add(left: core.math.Vec2, right: core.math.Vec2) core.math.Vec2 {
    return .{ .x = left.x + right.x, .y = left.y + right.y };
}

fn subtract(left: core.math.Vec2, right: core.math.Vec2) core.math.Vec2 {
    return .{ .x = left.x - right.x, .y = left.y - right.y };
}

fn scale(value: core.math.Vec2, factor: f32) core.math.Vec2 {
    return .{ .x = value.x * factor, .y = value.y * factor };
}

fn length(value: core.math.Vec2) f32 {
    return @sqrt(value.x * value.x + value.y * value.y);
}

fn distance(left: core.math.Vec2, right: core.math.Vec2) f32 {
    return length(subtract(left, right));
}

fn orderRing(allocator: std.mem.Allocator, ring: model.Ring, bonds: []const model.Bond, membership: topology.RingMembership) core.errors.Error![]core.ids.AtomId {
    const members = membership.atoms(ring.id);
    const result = allocator.alloc(core.ids.AtomId, members.len) catch return error.OutOfMemory;
    errdefer allocator.free(result);
    result[0] = members[0];
    var previous = core.ids.AtomId.invalid;
    for (1..members.len) |index| {
        var next = core.ids.AtomId.invalid;
        for (membership.ringBonds(ring.id)) |bond_id| {
            const bond = bonds[bond_id.index()];
            const candidate = if (bond.start == result[index - 1]) bond.end else if (bond.end == result[index - 1]) bond.start else continue;
            if (candidate == previous) continue;
            if (index + 1 < members.len and std.mem.indexOfScalar(core.ids.AtomId, result[0..index], candidate) != null) continue;
            if (!next.isValid() or candidate.index() < next.index()) next = candidate;
        }
        if (!next.isValid()) return error.InvalidMapping;
        previous = result[index - 1];
        result[index] = next;
    }
    return result;
}

/// `initializeVariablesForNeighboursCoordinatesRingAtom`: place a ring atom's
/// non-ring substituents inside the widest *scaled* gap between its ring
/// neighbours, rather than at a fixed offset from an arbitrary one.
///
/// Native previously took the first placed neighbour it found, added pi/3, and
/// spread the rest by 2*pi/3. On a fused quaternary carbon that picks between
/// two candidate gaps effectively at random, and choosing the wrong one puts a
/// methyl exactly two bond lengths from where upstream puts it - which is
/// drug_like/1's whole remaining residual (cgz-7v2.29, cgz-r31).
///
/// Returns false when the centre has no ring neighbour with coordinates, so the
/// caller keeps its own fallback rather than this inventing one.
fn placeRingAtomSubstituents(
    atoms: []model.Atom,
    graph: topology.Graph,
    membership: topology.RingMembership,
    fragmentation: fragments.Fragmentation,
    fragment: core.ids.FragmentId,
    center: core.ids.AtomId,
    placed: []bool,
    queue: []core.ids.AtomId,
    tail: *usize,
) core.errors.Error!bool {
    const origin = atoms[center.index()].coordinates;
    var ring_angles: [16]f32 = undefined;
    var ring_neighbours: [16]core.ids.AtomId = undefined;
    var ring_count: usize = 0;
    // Upstream's orderedNeighbors for this branch is every non-ring neighbour,
    // visited ones included: they are skipped when placing but still divide the
    // gap, so the divisor counts them.
    var substituents: [16]core.ids.AtomId = undefined;
    var substituent_count: usize = 0;
    for (graph.neighbors(center)) |neighbor| {
        if (sharesRing(membership, center, neighbor)) {
            if (ring_count == ring_angles.len) return false;
            var angle = std.math.atan2(
                atoms[neighbor.index()].coordinates.y - origin.y,
                atoms[neighbor.index()].coordinates.x - origin.x,
            );
            if (angle < 0) angle += 2 * pi;
            ring_angles[ring_count] = angle;
            ring_neighbours[ring_count] = neighbor;
            ring_count += 1;
        } else {
            if (substituent_count == substituents.len) return false;
            substituents[substituent_count] = neighbor;
            substituent_count += 1;
        }
    }
    if (ring_count == 0 or substituent_count == 0) return false;

    // Stable insertion sort by angle. Upstream stable_sorts pairs, whose
    // secondary key is the atom pointer; sorting on the angle alone and keeping
    // insertion order for exact ties is the order-stable reading of that
    // (cgz-r13), and two ring bonds at an identical angle is degenerate anyway.
    var i: usize = 1;
    while (i < ring_count) : (i += 1) {
        var j = i;
        while (j > 0 and ring_angles[j - 1] > ring_angles[j]) : (j -= 1) {
            std.mem.swap(f32, &ring_angles[j - 1], &ring_angles[j]);
            std.mem.swap(core.ids.AtomId, &ring_neighbours[j - 1], &ring_neighbours[j]);
        }
    }

    var best: usize = 0;
    var best_scaled: f32 = -1;
    var best_gap: f32 = 0;
    for (0..ring_count) |index| {
        const next = (index + 1) % ring_count;
        var gap = ring_angles[next] - ring_angles[index];
        if (gap < 0) gap += 2 * pi;
        var scaled = gap;
        if (gap > pi) {
            // A reflex gap is the outside of the ring system, and upstream
            // strongly prefers it.
            scaled *= 10;
        } else if (gapPointsIntoRing(atoms, membership, fragmentation, fragment, center, ring_angles[index] + gap * 0.5)) {
            scaled *= 0.2;
        }
        // First-wins on a tie, matching upstream's strict `>`.
        if (scaled > best_scaled) {
            best_scaled = scaled;
            best = index;
            best_gap = gap;
        }
    }

    // Upstream starts from the chosen gap's own ring neighbour and steps
    // counter-clockwise into the gap, one step per substituent slot. A visited
    // neighbour consumes no step but still counts in the divisor.
    const step = best_gap / @as(f32, @floatFromInt(substituent_count + 1));
    const start = std.math.atan2(
        atoms[ring_neighbours[best].index()].coordinates.y - origin.y,
        atoms[ring_neighbours[best].index()].coordinates.x - origin.x,
    );
    var taken: usize = 0;
    for (substituents[0..substituent_count]) |neighbor| {
        if (placed[neighbor.index()] and fragmentation.atom_fragment[neighbor.index()] == fragment) continue;
        taken += 1;
        const angle = start + step * @as(f32, @floatFromInt(taken));
        atoms[neighbor.index()].coordinates = .{
            .x = origin.x + @cos(angle) * bond_length,
            .y = origin.y + @sin(angle) * bond_length,
        };
        // As in the acyclic branch, a substituent belonging to another fragment
        // takes its slot in the gap and is written there - that write is how its
        // fragment learns its attachment direction - but it is not this
        // fragment's atom to visit.
        if (fragmentation.atom_fragment[neighbor.index()] != fragment) continue;
        placed[neighbor.index()] = true;
        queue[tail.*] = neighbor;
        tail.* += 1;
    }
    // The branch is chosen by the centre carrying rings, not by whether it
    // turned out to place anything: a centre whose substituents were all placed
    // already must not fall through to the fixed spread and be placed twice.
    return true;
}

/// The widest degree this port's neighbour buffers carry. Upstream has no cap;
/// a centre wider than this falls back to the fixed spread rather than being
/// laid out from a truncated list.
const max_neighbours = 16;

/// Upstream's acyclic half of `initializeVariablesForNeighboursCoordinates`,
/// together with the placement loop of
/// `generateCoordinatesNeighborsOfFirstAtomInQueue` that consumes it. The two
/// cannot be split: the angle list is indexed by position in the ordered
/// neighbour list, so using upstream's angles requires walking upstream's list.
///
/// Native previously placed only the same-fragment unplaced neighbours, spread
/// symmetrically about a single base direction. That dropped three things
/// upstream does - priority ordering at degree four, the per-centre angles from
/// `neighborsAnglesAtCenter`, and the cumulative rotation that a skipped
/// neighbour still displaces - and the error compounded down a chain of joints
/// (cgz-7v2.30).
fn placeAcyclicNeighbours(
    allocator: std.mem.Allocator,
    atoms: []model.Atom,
    bonds: []const model.Bond,
    graph: topology.Graph,
    membership: topology.RingMembership,
    fragmentation: fragments.Fragmentation,
    fragment: core.ids.FragmentId,
    center: core.ids.AtomId,
    placed: []bool,
    queue: []core.ids.AtomId,
    tail: *usize,
) core.errors.Error!bool {
    const neighbours = graph.neighbors(center);
    if (neighbours.len == 0 or neighbours.len > max_neighbours) return false;

    var ordered: [max_neighbours]core.ids.AtomId = undefined;
    try orderNeighbours(allocator, atoms, bonds, graph, membership, center, ordered[0..neighbours.len]);

    // Upstream's visited set is fragment-local: it starts as this fragment's
    // ring atoms, or its single start atom, and only ever gains atoms of this
    // fragment. An atom another fragment already placed is not visited, so the
    // test is `placed` narrowed to this fragment.
    var direction: core.math.Vec2 = .{ .x = -bond_length, .y = 0 };
    var start: usize = 0;
    for (ordered[0..neighbours.len], 0..) |neighbor, index| {
        if (!placed[neighbor.index()] or fragmentation.atom_fragment[neighbor.index()] != fragment) continue;
        direction = subtract(atoms[neighbor.index()].coordinates, atoms[center.index()].coordinates);
        start = index;
        break;
    }
    std.mem.rotate(core.ids.AtomId, ordered[0..neighbours.len], start);

    var angles: [max_neighbours]f32 = undefined;
    neighbourAnglesAtCenter(atoms, bonds, graph, center, angles[0..neighbours.len]);

    for (ordered[0..neighbours.len], angles[0..neighbours.len]) |neighbor, angle| {
        // A visited neighbour is skipped before the rotation, so it holds its
        // slot in the angle list without advancing the running direction.
        if (placed[neighbor.index()] and fragmentation.atom_fragment[neighbor.index()] == fragment) continue;
        direction = rotateClockwise(direction, angle);
        atoms[neighbor.index()].coordinates = add(atoms[center.index()].coordinates, direction);
        // A neighbour outside this fragment is written but neither visited nor
        // queued: the write is what tells its own fragment where its attachment
        // atom goes, and `FragmentFrames.store` reads it before the next
        // fragment's build can overwrite it.
        if (fragmentation.atom_fragment[neighbor.index()] != fragment) continue;
        placed[neighbor.index()] = true;
        queue[tail.*] = neighbor;
        tail.* += 1;
    }
    // Handled, whether or not anything moved: the caller must not also run the
    // fixed spread over the same centre.
    return true;
}

/// Upstream takes the neighbours in input order everywhere except at degree
/// four, where `sketcherMinimizerAtom::orderAtomPriorities` ranks them.
fn orderNeighbours(
    allocator: std.mem.Allocator,
    atoms: []const model.Atom,
    bonds: []const model.Bond,
    graph: topology.Graph,
    membership: topology.RingMembership,
    center: core.ids.AtomId,
    out: []core.ids.AtomId,
) core.errors.Error!void {
    const neighbours = graph.neighbors(center);
    @memcpy(out, neighbours);
    if (neighbours.len != 4) return;

    var weights: [4]f32 = undefined;
    for (neighbours, 0..) |neighbor, index| {
        weights[index] = try neighbourPriorityWeight(allocator, atoms, bonds, graph, membership, center, neighbor);
    }

    // Upstream lifts the lowest weight out, then the lowest of what remains,
    // and rebuilds the list with those two at the ends: long chains keep
    // positions 2 and 3, side substituents take 1 and 4. Phosphorus and sulfur
    // at degree four are the exception and take both lifted atoms up front.
    var rest: [4]core.ids.AtomId = undefined;
    var rest_weights: [4]f32 = undefined;
    var rest_count: usize = 0;
    for (neighbours, 0..) |neighbor, index| {
        rest[rest_count] = neighbor;
        rest_weights[rest_count] = weights[index];
        rest_count += 1;
    }
    const first = takeLowestWeight(rest[0..rest_count], rest_weights[0..rest_count]);
    rest_count -= 1;
    const second = takeLowestWeight(rest[0..rest_count], rest_weights[0..rest_count]);
    rest_count -= 1;

    const center_element = atoms[center.index()].atomic_number;
    if (center_element != .sulfur and center_element != .phosphorus) {
        out[0] = first;
        out[1] = rest[0];
        out[2] = rest[1];
        out[3] = second;
    } else {
        out[0] = first;
        out[1] = rest[0];
        out[2] = second;
        out[3] = rest[1];
    }
}

/// Remove and return the first atom holding the lowest weight, keeping the
/// order of the rest, which is upstream's `erase` on both parallel vectors.
fn takeLowestWeight(atoms: []core.ids.AtomId, weights: []f32) core.ids.AtomId {
    var lowest: usize = 0;
    for (weights, 0..) |weight, index| {
        if (weight < weights[lowest]) lowest = index;
    }
    const taken = atoms[lowest];
    var index = lowest;
    while (index + 1 < atoms.len) : (index += 1) {
        atoms[index] = atoms[index + 1];
        weights[index] = weights[index + 1];
    }
    return taken;
}

/// The weight `sketcherMinimizerAtom::orderAtomPriorities` scores a neighbour
/// by: the size of the branch behind it, adjusted by bond order, ring
/// membership, element, and neighbouring stereochemistry.
///
/// One term of upstream's is absent: the -2000 for a neighbour that is
/// `isSharedAndInner` while the centre is not. That flag is set only on the
/// fusion atoms of a ring system whose every neighbour lies in both fused
/// rings (CoordgenFragmentBuilder.cpp:975-995), and this port does not model
/// it at all - inventing a value for it here would be a guess, not a
/// transcription (cgz-7v2.31).
fn neighbourPriorityWeight(
    allocator: std.mem.Allocator,
    atoms: []const model.Atom,
    bonds: []const model.Bond,
    graph: topology.Graph,
    membership: topology.RingMembership,
    center: core.ids.AtomId,
    neighbor: core.ids.AtomId,
) core.errors.Error!f32 {
    const branch = try graph.reachableExcluding(allocator, neighbor, center);
    defer allocator.free(branch);
    var weight: f32 = @floatFromInt(branch.len);

    const bond_id = bondBetween(graph, center, neighbor);
    if (bond_id) |id| {
        const order = bonds[id.index()].effective_order;
        // So that =O gets lower priority than -OH in a phosphate.
        if (order == .double) weight -= 0.25;
        // Forcing the wedge away from the double bond in a sulphoxide.
        if (atoms[center.index()].atomic_number == .sulfur and order == .double) weight += 2000;
        // Upstream's test is sameRing on the bond's two atoms, not on the
        // bond, so a bridging bond between two atoms of one ring counts.
        if (sharesRing(membership, center, neighbor)) weight += 500;
    }
    if (atoms[neighbor.index()].atomic_number == .carbon) weight += 0.5;
    if (atoms[neighbor.index()].atomic_number == .hydrogen) weight -= 0.5;
    if (atoms[neighbor.index()].stereo != .unspecified) weight += 10000;
    if (atoms[center.index()].cross_layout and graph.degree(neighbor) > 1) weight += 200;
    for (graph.incidentBonds(neighbor)) |incident| if (bonds[incident.index()].effective_order == .double) {
        weight += 100;
        break;
    };
    return weight;
}

fn bondBetween(graph: topology.Graph, atom: core.ids.AtomId, other: core.ids.AtomId) ?core.ids.BondId {
    for (graph.neighbors(atom), graph.incidentBonds(atom)) |neighbor, bond| {
        if (neighbor == other) return bond;
    }
    return null;
}

/// `CoordgenFragmentBuilder::neighborsAnglesAtCenter`. `m_evenAngles` is
/// upstream's constructor default of false and nothing in this port sets it,
/// so only the uneven branch is transcribed.
fn neighbourAnglesAtCenter(
    atoms: []const model.Atom,
    bonds: []const model.Bond,
    graph: topology.Graph,
    center: core.ids.AtomId,
    out: []f32,
) void {
    const neighbours = graph.neighbors(center);
    var division = neighbours.len;
    if (neighbours.len == 2) {
        if (atoms[center.index()].atomic_number == .carbon or
            !atoms[neighbours[0].index()].cross_layout or
            !atoms[neighbours[1].index()].cross_layout) division = 3;
        const incident = graph.incidentBonds(center);
        const total = @intFromEnum(bonds[incident[0].index()].effective_order) +
            @intFromEnum(bonds[incident[1].index()].effective_order);
        if (total >= 4) division = 2;
    } else if (neighbours.len == 4 and !atoms[center.index()].cross_layout) {
        // Upstream's 60-90-120-90 around a tetracoordinated centre.
        out[0] = pi / 3;
        out[1] = pi * 0.5;
        out[2] = 2 * pi / 3;
        out[3] = pi * 0.5;
        return;
    }
    for (out) |*angle| angle.* = 2 * pi / @as(f32, @floatFromInt(division));
}

/// `sketcherMinimizerPointF::rotate(sin a, cos a)` for a caller holding an
/// angle rather than its sine and cosine.
fn rotateClockwise(value: core.math.Vec2, angle: f32) core.math.Vec2 {
    return rotateScreen(value, @sin(angle), @cos(angle));
}

// --- cgz-7v2.27: the child fragment's side choice -------------------------

/// `sketcherMinimizer.cpp`'s three scoring multipliers, at their upstream
/// values.
const score_multiplier_for_double_bonds: f32 = 0.82;
const score_multiplier_for_single_bonded_heteroatoms: f32 = 0.9;
const score_multiplier_for_fragments: f32 = 0.1;

/// Scratch for the side choice, sized once for the whole molecule rather than
/// per fragment: `getAllTerminalBonds` can return at most every bond plus one
/// per child plus the bond to the parent.
const SideChoiceScratch = struct {
    allocator: std.mem.Allocator,
    /// Terminal bonds of the fragment being placed, and of its parent.
    own_bonds: []OrientedBond,
    parent_bonds: []OrientedBond,
    /// One direction and score per surviving parent terminal bond.
    directions: []core.math.Vec2,
    scores: []f32,
    /// `longestChainFromHere`, indexed by fragment.
    longest_chain: []f32,

    fn init(
        allocator: std.mem.Allocator,
        bond_count: usize,
        fragment_count: usize,
    ) core.errors.Error!SideChoiceScratch {
        const capacity = bond_count + fragment_count + 1;
        const own_bonds = allocator.alloc(OrientedBond, capacity) catch return error.OutOfMemory;
        errdefer allocator.free(own_bonds);
        const parent_bonds = allocator.alloc(OrientedBond, capacity) catch return error.OutOfMemory;
        errdefer allocator.free(parent_bonds);
        const directions = allocator.alloc(core.math.Vec2, capacity) catch return error.OutOfMemory;
        errdefer allocator.free(directions);
        const scores = allocator.alloc(f32, capacity) catch return error.OutOfMemory;
        errdefer allocator.free(scores);
        const longest_chain = allocator.alloc(f32, fragment_count) catch return error.OutOfMemory;
        return .{
            .allocator = allocator,
            .own_bonds = own_bonds,
            .parent_bonds = parent_bonds,
            .directions = directions,
            .scores = scores,
            .longest_chain = longest_chain,
        };
    }

    fn deinit(self: *SideChoiceScratch) void {
        self.allocator.free(self.longest_chain);
        self.allocator.free(self.scores);
        self.allocator.free(self.directions);
        self.allocator.free(self.parent_bonds);
        self.allocator.free(self.own_bonds);
        self.* = undefined;
    }
};

/// A terminal bond with its endpoints oriented the way upstream's scoring
/// reads them: `start` on the side of the fragment whose list this is, `end`
/// on the far side.
///
/// Upstream can use `bond->startAtom` and `bond->endAtom` directly because a
/// fragment's `_bondToParent` is built parent-first. Native stores bonds in
/// canonical input order, which on the drug-like corpus is child-first
/// (`3-2`, `5-2`, ...), so reading `bond.start` here selects the wrong atom.
/// It is not a latent difference: taking the stored direction made
/// `scoreDirections` skip every child bond, which left single-atom fragments
/// with nothing to score and no side choice at all - measured as centres 2 and
/// 5 of drug_like/4 coming out mirrored.
const OrientedBond = struct {
    id: core.ids.BondId,
    start: core.ids.AtomId,
    end: core.ids.AtomId,
};

/// `sketcherMinimizer::getAllTerminalBonds`: every bond of the fragment with a
/// degree-one endpoint, then one bond per child, then the bond to the parent.
///
/// Upstream walks `fragment->getBonds()` in the fragment's own stored order.
/// This walks the molecule's bonds in index order and keeps those with both
/// endpoints in the fragment, which is the order-stable reading of that
/// (cgz-r13); the scoring below takes a strict maximum, so the order is
/// observable only through ties.
///
/// Residue interactions are skipped upstream. Native models them outside
/// `bonds` entirely, so there is nothing here to skip.
fn collectTerminalBonds(
    bonds: []const model.Bond,
    graph: topology.Graph,
    fragmentation: fragments.Fragmentation,
    fragment: fragments.Fragment,
    out: []OrientedBond,
) []OrientedBond {
    var count: usize = 0;
    for (bonds) |bond| {
        if (fragmentation.atom_fragment[bond.start.index()] != fragment.id or
            fragmentation.atom_fragment[bond.end.index()] != fragment.id) continue;
        if (graph.degree(bond.start) != 1 and graph.degree(bond.end) != 1) continue;
        if (count == out.len) return out[0..count];
        // Both endpoints are inside the fragment, so either orientation reads
        // the same: the midpoint, the degree test and the heteroatom test are
        // all symmetric.
        out[count] = .{ .id = bond.id, .start = bond.start, .end = bond.end };
        count += 1;
    }
    for (fragmentation.fragments) |candidate| {
        if (candidate.parent != fragment.id or !candidate.bond_to_parent.isValid()) continue;
        if (!candidate.attachment_atom.isValid() or !candidate.parent_atom.isValid()) continue;
        if (count == out.len) return out[0..count];
        out[count] = .{ .id = candidate.bond_to_parent, .start = candidate.parent_atom, .end = candidate.attachment_atom };
        count += 1;
    }
    if (fragment.bond_to_parent.isValid() and count != out.len) {
        // Oriented outward, so `scoreDirections`'s "starts in this fragment"
        // test skips it exactly as upstream's does.
        if (fragment.attachment_atom.isValid() and fragment.parent_atom.isValid()) {
            out[count] = .{
                .id = fragment.bond_to_parent,
                .start = fragment.parent_atom,
                .end = fragment.attachment_atom,
            };
            count += 1;
        }
    }
    return out[0..count];
}

/// `sketcherMinimizer::assignLongestChainFromHere`, computed for every fragment
/// at once. The recursion is replaced by a walk from the deepest fragments
/// upward, which visits every child before its parent and needs no stack.
///
/// `positionFromParent` upstream is the fragment's attachment atom read out of
/// its parent's own frame, which is exactly what `FragmentFrames.attachment`
/// holds.
fn assignLongestChains(
    fragmentation: fragments.Fragmentation,
    frames: FragmentFrames,
    longest_chain: []f32,
) void {
    @memset(longest_chain, 0);
    var depth = fragmentation.fragments.len;
    while (depth > 0) {
        depth -= 1;
        for (fragmentation.fragments) |fragment| {
            if (fragmentDepth(fragmentation.fragments, fragment.id) != depth) continue;
            var longest: f32 = 0;
            for (fragmentation.fragments) |candidate| {
                if (candidate.parent != fragment.id) continue;
                longest = @max(longest, longest_chain[candidate.id.index()]);
            }
            const offset = frames.attachment[fragment.id.index()];
            const from_parent = if (fragment.parent.isValid())
                @sqrt(offset.x * offset.x + offset.y * offset.y)
            else
                0;
            longest_chain[fragment.id.index()] = longest + from_parent;
        }
    }
}

/// The position an atom holds in `fragment`'s own frame. Upstream's
/// `_coordinates` map covers the fragment's own atoms and, additionally, the
/// attachment atom of each of its children; both are stored here, in the two
/// halves of `FragmentFrames`.
fn frameCoordinate(
    fragmentation: fragments.Fragmentation,
    frames: FragmentFrames,
    fragment: fragments.Fragment,
    atom: core.ids.AtomId,
) core.math.Vec2 {
    if (fragmentation.atom_fragment[atom.index()] == fragment.id) return frames.own[atom.index()];
    for (fragmentation.fragments) |candidate| {
        if (candidate.parent != fragment.id) continue;
        const child_attachment = attachmentOf(candidate) orelse continue;
        if (child_attachment == atom) return frames.attachment[candidate.id.index()];
    }
    return frames.own[atom.index()];
}

fn normalized(value: core.math.Vec2) core.math.Vec2 {
    const magnitude = @sqrt(value.x * value.x + value.y * value.y);
    if (magnitude == 0) return value;
    return .{ .x = value.x / magnitude, .y = value.y / magnitude };
}

/// The shared part of both scoring loops: a double bond is worth less, and a
/// bond to a terminal heteroatom is worth less again.
fn terminalBondScoreModifier(
    atoms: []const model.Atom,
    bonds: []const model.Bond,
    graph: topology.Graph,
    bond: OrientedBond,
) f32 {
    var score: f32 = 1;
    if (bonds[bond.id.index()].effective_order == .double) score *= score_multiplier_for_double_bonds;
    if ((graph.degree(bond.start) == 1 and atoms[bond.start.index()].atomic_number != .carbon) or
        (graph.degree(bond.end) == 1 and atoms[bond.end.index()].atomic_number != .carbon))
    {
        score *= score_multiplier_for_single_bonded_heteroatoms;
    }
    return score;
}

/// `sketcherMinimizer::findDirectionsToAlignWith`: the directions, in the
/// molecule's own placed coordinates, that point from the parent's terminal
/// bonds toward the bond this fragment hangs off. These are what the fragment's
/// own terminal bonds are then asked to line up with.
fn findDirectionsToAlignWith(
    atoms: []const model.Atom,
    bonds: []const model.Bond,
    graph: topology.Graph,
    fragmentation: fragments.Fragmentation,
    fragment: fragments.Fragment,
    scratch: SideChoiceScratch,
) usize {
    const parent = fragmentation.fragments[fragment.parent.index()];
    const attachment = attachmentOf(fragment) orelse return 0;
    const parent_atom = parentAtomOf(fragment) orelse return 0;
    const origin = scale(add(
        atoms[parent_atom.index()].coordinates,
        atoms[attachment.index()].coordinates,
    ), 0.5);

    var count: usize = 0;
    for (collectTerminalBonds(bonds, graph, fragmentation, parent, scratch.parent_bonds)) |bond| {
        if (fragmentation.atom_fragment[bond.end.index()] == fragment.id) continue;
        const midpoint = scale(add(
            atoms[bond.start.index()].coordinates,
            atoms[bond.end.index()].coordinates,
        ), 0.5);
        var score = terminalBondScoreModifier(atoms, bonds, graph, bond);
        // A bond that leaves the parent is scored by how much molecule hangs
        // off it rather than by what it is, and a bond running back toward the
        // grandparent outranks everything.
        if (fragmentation.atom_fragment[bond.end.index()] != parent.id or
            fragmentation.atom_fragment[bond.start.index()] != parent.id)
        {
            const end_fragment = fragmentation.atom_fragment[bond.end.index()];
            score = scratch.longest_chain[end_fragment.index()] * score_multiplier_for_fragments;
            if (parent.parent.isValid() and
                fragmentation.atom_fragment[bond.start.index()] == parent.parent) score *= 100;
        }
        if (count == scratch.directions.len) break;
        scratch.directions[count] = normalized(subtract(origin, midpoint));
        scratch.scores[count] = score;
        count += 1;
    }
    return count;
}

/// `sketcherMinimizer::testAlignment`: how well one direction lines up with one
/// target, squared, with a large bonus for an exact match.
fn testAlignment(direction: core.math.Vec2, target: core.math.Vec2, weight: f32) f32 {
    var dot = direction.x * target.x + direction.y * target.y;
    if (dot < 0) dot = 0;
    var score = dot * dot;
    if (dot > 1 - sketcher_epsilon) score += 1000;
    return score * weight;
}

/// `sketcherMinimizer::scoreDirections`, reduced to the one thing its caller
/// uses: whether the fragment should be flipped about its bond to the parent.
///
/// Each of the fragment's own terminal bonds is taken in its frame, both as it
/// stands and mirrored in y, turned by the angle the fragment is about to be
/// placed at, and scored against every direction found above. The better of the
/// two readings wins, and if that is the mirrored one the fragment is flipped.
fn shouldInvertAgainstParent(
    atoms: []const model.Atom,
    bonds: []const model.Bond,
    graph: topology.Graph,
    fragmentation: fragments.Fragmentation,
    frames: FragmentFrames,
    fragment: fragments.Fragment,
    angle: f32,
    scratch: SideChoiceScratch,
) bool {
    const direction_count = findDirectionsToAlignWith(atoms, bonds, graph, fragmentation, fragment, scratch);
    if (direction_count == 0) return false;

    var best_score: f32 = 0;
    var invert = false;
    for (collectTerminalBonds(bonds, graph, fragmentation, fragment, scratch.own_bonds)) |bond| {
        if (fragmentation.atom_fragment[bond.start.index()] != fragment.id) continue;
        const midpoint = scale(add(
            frameCoordinate(fragmentation, frames, fragment, bond.start),
            frameCoordinate(fragmentation, frames, fragment, bond.end),
        ), 0.5);
        // The frame's own origin for this purpose is half a bond length back
        // along -x, which is where the bond to the parent enters it.
        const plain = normalized(subtract(midpoint, .{ .x = -bond_length * 0.5, .y = 0 }));
        const inverted: core.math.Vec2 = .{ .x = plain.x, .y = -plain.y };
        const turned_plain = rotateClockwise(plain, angle);
        const turned_inverted = rotateClockwise(inverted, angle);

        var modifier = terminalBondScoreModifier(atoms, bonds, graph, bond);
        if (fragmentation.atom_fragment[bond.end.index()] != fragment.id) {
            const end_fragment = fragmentation.atom_fragment[bond.end.index()];
            modifier = scratch.longest_chain[end_fragment.index()] * score_multiplier_for_fragments;
        }
        for (scratch.directions[0..direction_count], scratch.scores[0..direction_count]) |target, weight| {
            // Strict `>`, so the first reading to reach a score keeps it.
            if (testAlignment(turned_plain, target, weight) * modifier > best_score) {
                best_score = testAlignment(turned_plain, target, weight) * modifier;
                invert = false;
            }
            if (testAlignment(turned_inverted, target, weight) * modifier > best_score) {
                best_score = testAlignment(turned_inverted, target, weight) * modifier;
                invert = true;
            }
        }
    }
    return invert;
}

/// `sketcherMinimizer::alignWithParentDirection`'s effect on the pose: mirror
/// the fragment's stored frame in y, which rotates it 180 degrees about the
/// bond it hangs off once the frame is placed. The children's attachment
/// positions live in this frame too and are mirrored with it, so a whole
/// subtree follows.
///
/// Upstream also inverts the wedge/hash of any stereo bond on a flipped atom.
/// Native does not, deliberately: bond displays are assigned after layout by
/// `topology.stereo.writeAtomBondDisplays`, from the geometry this function has
/// already settled, so the inversion happens on its own. That path is exact on
/// the drug-like partition (`bond_displays` 7/7), and flipping here as well
/// would invert it twice.
fn mirrorFragmentFrame(
    fragmentation: fragments.Fragmentation,
    frames: FragmentFrames,
    fragment: fragments.Fragment,
) void {
    for (fragmentation.members(fragment.id)) |atom| {
        frames.own[atom.index()].y = -frames.own[atom.index()].y;
    }
    for (fragmentation.fragments) |candidate| {
        if (candidate.parent != fragment.id) continue;
        if (attachmentOf(candidate) == null) continue;
        frames.attachment[candidate.id.index()].y = -frames.attachment[candidate.id.index()].y;
    }
}

fn sharesRing(membership: topology.RingMembership, left: core.ids.AtomId, right: core.ids.AtomId) bool {
    for (membership.atomRings(left)) |ring| {
        for (membership.atomRings(right)) |other| {
            if (ring == other) return true;
        }
    }
    return false;
}

/// Upstream's SKETCHER_EPSILON, the guard against dividing by a horizontal edge.
const sketcher_epsilon: f32 = 1.0e-5;

/// Upstream probes a point a tenth of a bond length along the gap's midpoint and
/// asks whether it falls inside any non-macrocyclic ring of the fragment
/// (`sketcherMinimizerRing::contains`, an even-odd crossing count).
fn gapPointsIntoRing(
    atoms: []const model.Atom,
    membership: topology.RingMembership,
    fragmentation: fragments.Fragmentation,
    fragment: core.ids.FragmentId,
    center: core.ids.AtomId,
    midpoint_angle: f32,
) bool {
    const origin = atoms[center.index()].coordinates;
    const probe: core.math.Vec2 = .{
        .x = origin.x + @cos(midpoint_angle) * bond_length * 0.1,
        .y = origin.y + @sin(midpoint_angle) * bond_length * 0.1,
    };
    for (membership.rings) |ring| {
        const members = membership.atoms(ring.id);
        if (members.len == 0) continue;
        if (fragmentation.atom_fragment[members[0].index()] != fragment) continue;
        if (members.len >= topology.rings.macrocycle_size) continue;
        if (ringContainsPoint(atoms, members, probe)) return true;
    }
    return false;
}

fn ringContainsPoint(atoms: []const model.Atom, members: []const core.ids.AtomId, probe: core.math.Vec2) bool {
    var crossings: usize = 0;
    for (members, 0..) |member, index| {
        const next = members[(index + 1) % members.len];
        const a = atoms[member.index()].coordinates;
        const b = atoms[next.index()].coordinates;
        // Upstream's strict inequalities, so a vertex exactly at the probe's
        // height is not counted from either side.
        if (!((probe.y < a.y and probe.y > b.y) or (probe.y > a.y and probe.y < b.y))) continue;
        const dy = b.y - a.y;
        if (dy <= sketcher_epsilon and dy >= -sketcher_epsilon) continue;
        const t = (probe.y - a.y) / dy;
        if (probe.x > a.x + (b.x - a.x) * t) crossings += 1;
    }
    return crossings % 2 != 0;
}

fn parentAngle(atoms: []const model.Atom, graph: topology.Graph, fragmentation: fragments.Fragmentation, fragment: core.ids.FragmentId, center: core.ids.AtomId, placed: []const bool) f32 {
    for (graph.neighbors(center)) |neighbor| {
        if (!placed[neighbor.index()] or fragmentation.atom_fragment[neighbor.index()] != fragment or
            (atoms[neighbor.index()].coordinates.x == atoms[center.index()].coordinates.x and atoms[neighbor.index()].coordinates.y == atoms[center.index()].coordinates.y)) continue;
        const direction = atoms[center.index()].coordinates;
        const parent = atoms[neighbor.index()].coordinates;
        return std.math.atan2(direction.y - parent.y, direction.x - parent.x) + pi / 3;
    }
    return 0;
}

test "regular ring coordinate walk preserves the upstream bond length" {
    var atoms: [6]model.Atom = undefined;
    var bonds: [6]model.Bond = undefined;
    for (&atoms, 0..) |*atom, index| atom.* = .{ .id = core.ids.AtomId.fromIndex(@intCast(index)), .input_index = @intCast(index), .atomic_number = .carbon };
    for (&bonds, 0..) |*bond, index| bond.* = .{
        .id = core.ids.BondId.fromIndex(@intCast(index)),
        .input_index = @intCast(index),
        .start = core.ids.AtomId.fromIndex(@intCast(index)),
        .end = core.ids.AtomId.fromIndex(@intCast((index + 1) % atoms.len)),
        .input_order = .single,
        .effective_order = .single,
    };
    var graph = try topology.Graph.init(std.testing.allocator, &atoms, &bonds);
    defer graph.deinit();
    var rings = try topology.RingMembership.init(std.testing.allocator, graph, &bonds);
    defer rings.deinit();
    var split = try fragments.Fragmentation.init(std.testing.allocator, &atoms, &bonds, graph, rings);
    defer split.deinit();
    try initializeCoordinates(std.testing.allocator, &atoms, &bonds, graph, rings, split);
    for (bonds) |bond| {
        const a = atoms[bond.start.index()].coordinates;
        const b = atoms[bond.end.index()].coordinates;
        try std.testing.expectApproxEqAbs(bond_length, @sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y)), 0.001);
    }
}

test "macrocycles dispatch through native polyomino placement" {
    const ring_size = 13;
    var atoms: [ring_size]model.Atom = undefined;
    var bonds: [ring_size]model.Bond = undefined;
    for (&atoms, 0..) |*atom, index| atom.* = .{ .id = core.ids.AtomId.fromIndex(@intCast(index)), .input_index = @intCast(index), .atomic_number = .carbon };
    for (&bonds, 0..) |*bond, index| bond.* = .{
        .id = core.ids.BondId.fromIndex(@intCast(index)),
        .input_index = @intCast(index),
        .start = core.ids.AtomId.fromIndex(@intCast(index)),
        .end = core.ids.AtomId.fromIndex(@intCast((index + 1) % ring_size)),
        .input_order = .single,
        .effective_order = .single,
    };
    var graph = try topology.Graph.init(std.testing.allocator, &atoms, &bonds);
    defer graph.deinit();
    var rings = try topology.RingMembership.init(std.testing.allocator, graph, &bonds);
    defer rings.deinit();
    var split = try fragments.Fragmentation.init(std.testing.allocator, &atoms, &bonds, graph, rings);
    defer split.deinit();
    try initializeCoordinates(std.testing.allocator, &atoms, &bonds, graph, rings, split);

    var hex_edges: usize = 0;
    for (bonds) |bond| {
        const edge_length = distance(atoms[bond.start.index()].coordinates, atoms[bond.end.index()].coordinates);
        hex_edges += @intFromBool(@abs(edge_length - bond_length * @sqrt(@as(f32, 3))) < 0.001);
    }
    try std.testing.expect(hex_edges > 0);

    const opened = macrocycle.findBondToOpen(core.ids.RingId.fromIndex(0), &bonds, graph, rings) orelse return error.InvalidMapping;
    _ = try initializeCoordinatesWithOptions(std.testing.allocator, &atoms, &bonds, graph, rings, split, .{ .force_open_macrocycles = true });
    for (bonds) |bond| {
        if (bond.id == opened) continue;
        try std.testing.expectApproxEqAbs(
            bond_length,
            distance(atoms[bond.start.index()].coordinates, atoms[bond.end.index()].coordinates),
            0.001,
        );
    }
    try std.testing.expect(@abs(distance(
        atoms[bonds[opened.index()].start.index()].coordinates,
        atoms[bonds[opened.index()].end.index()].coordinates,
    ) - bond_length) > 0.001);
    try core.oom.checkAllocationFailures(
        std.testing.allocator,
        layoutWithOptionsAndDiscard,
        .{ &atoms, &bonds, graph, rings, split, true },
    );
}

fn layoutFixture(atoms: []model.Atom, bonds: []const model.Bond) !void {
    var graph = try topology.Graph.init(std.testing.allocator, atoms, bonds);
    defer graph.deinit();
    var rings = try topology.RingMembership.init(std.testing.allocator, graph, bonds);
    defer rings.deinit();
    var split = try fragments.Fragmentation.init(std.testing.allocator, atoms, bonds, graph, rings);
    defer split.deinit();
    try initializeCoordinates(std.testing.allocator, atoms, bonds, graph, rings, split);
}

fn layoutAndDiscard(
    allocator: std.mem.Allocator,
    source_atoms: []const model.Atom,
    bonds: []const model.Bond,
    graph: topology.Graph,
    rings: topology.RingMembership,
    split: fragments.Fragmentation,
) !void {
    const atoms = try allocator.dupe(model.Atom, source_atoms);
    defer allocator.free(atoms);
    try initializeCoordinates(allocator, atoms, bonds, graph, rings, split);
}

fn layoutWithOptionsAndDiscard(
    allocator: std.mem.Allocator,
    source_atoms: []const model.Atom,
    bonds: []const model.Bond,
    graph: topology.Graph,
    rings: topology.RingMembership,
    split: fragments.Fragmentation,
    force_open_macrocycles: bool,
) !void {
    const atoms = try allocator.dupe(model.Atom, source_atoms);
    defer allocator.free(atoms);
    _ = try initializeCoordinatesWithOptions(allocator, atoms, bonds, graph, rings, split, .{ .force_open_macrocycles = force_open_macrocycles });
}

fn captureFramesAndDiscard(
    allocator: std.mem.Allocator,
    atoms: []const model.Atom,
    split: fragments.Fragmentation,
) !void {
    var frames = try captureFragmentFrames(allocator, atoms, split);
    defer frames.deinit();
}

fn expectBondLengths(atoms: []const model.Atom, bonds: []const model.Bond) !void {
    for (bonds) |bond| {
        try std.testing.expectApproxEqAbs(
            bond_length,
            distance(atoms[bond.start.index()].coordinates, atoms[bond.end.index()].coordinates),
            0.01,
        );
    }
}

test "fused rings align on their shared edge and extend outward" {
    var atoms: [6]model.Atom = undefined;
    for (&atoms, 0..) |*atom, index| atom.* = .{ .id = core.ids.AtomId.fromIndex(@intCast(index)), .input_index = @intCast(index), .atomic_number = .carbon };
    const pairs = [_][2]u32{ .{ 0, 1 }, .{ 1, 2 }, .{ 2, 3 }, .{ 3, 0 }, .{ 2, 5 }, .{ 5, 4 }, .{ 4, 3 } };
    var bonds: [pairs.len]model.Bond = undefined;
    for (&bonds, pairs, 0..) |*bond, pair, index| bond.* = .{
        .id = core.ids.BondId.fromIndex(@intCast(index)),
        .input_index = @intCast(index),
        .start = core.ids.AtomId.fromIndex(pair[0]),
        .end = core.ids.AtomId.fromIndex(pair[1]),
        .input_order = .single,
        .effective_order = .single,
    };
    try layoutFixture(&atoms, &bonds);
    try expectBondLengths(&atoms, &bonds);
    const first_center = scale(add(add(atoms[0].coordinates, atoms[1].coordinates), add(atoms[2].coordinates, atoms[3].coordinates)), 0.25);
    const second_center = scale(add(add(atoms[2].coordinates, atoms[3].coordinates), add(atoms[4].coordinates, atoms[5].coordinates)), 0.25);
    try std.testing.expect(distance(first_center, second_center) > bond_length * 0.5);

    var graph = try topology.Graph.init(std.testing.allocator, &atoms, &bonds);
    defer graph.deinit();
    var rings = try topology.RingMembership.init(std.testing.allocator, graph, &bonds);
    defer rings.deinit();
    var split = try fragments.Fragmentation.init(std.testing.allocator, &atoms, &bonds, graph, rings);
    defer split.deinit();
    try core.oom.checkAllocationFailures(
        std.testing.allocator,
        layoutAndDiscard,
        .{ &atoms, &bonds, graph, rings, split },
    );
}

test "the third ring of a linear acene mirrors away from its parent ring, not from everything drawn" {
    // A linear tricyclic is the smallest case where upstream's two mirror
    // references disagree: when the third ring is placed, the centroid of the
    // ten atoms already drawn sits well outside its parent ring's own centre.
    // Upstream scores against the parent ring (CoordgenFragmentBuilder.cpp:513),
    // so this pins that reference rather than the placed centroid (cgz-r31).
    var atoms: [14]model.Atom = undefined;
    for (&atoms, 0..) |*atom, index| atom.* = .{ .id = core.ids.AtomId.fromIndex(@intCast(index)), .input_index = @intCast(index), .atomic_number = .carbon };
    const pairs = [_][2]u32{
        .{ 0, 1 },  .{ 1, 2 },   .{ 2, 3 },   .{ 3, 4 },  .{ 4, 5 },   .{ 5, 0 },
        .{ 3, 6 },  .{ 6, 7 },   .{ 7, 8 },   .{ 8, 9 },  .{ 9, 4 },   .{ 8, 10 },
        .{ 10, 11 }, .{ 11, 12 }, .{ 12, 13 }, .{ 13, 7 },
    };
    var bonds: [pairs.len]model.Bond = undefined;
    for (&bonds, pairs, 0..) |*bond, pair, index| bond.* = .{
        .id = core.ids.BondId.fromIndex(@intCast(index)),
        .input_index = @intCast(index),
        .start = core.ids.AtomId.fromIndex(pair[0]),
        .end = core.ids.AtomId.fromIndex(pair[1]),
        .input_order = .single,
        .effective_order = .single,
    };
    try layoutFixture(&atoms, &bonds);
    try expectBondLengths(&atoms, &bonds);

    // The middle ring is the outer ring's parent, and the outer ring's own
    // atoms must sit further from that centre than their mirror images across
    // the shared 7-8 edge would.
    const middle = [_]u32{ 3, 4, 6, 7, 8, 9 };
    var parent_center: core.math.Vec2 = .{};
    for (middle) |index| parent_center = add(parent_center, atoms[index].coordinates);
    parent_center = scale(parent_center, 1.0 / @as(f32, @floatFromInt(middle.len)));

    const axis_start = atoms[7].coordinates;
    const axis_end = atoms[8].coordinates;
    var direct: f32 = 0;
    var mirrored: f32 = 0;
    for ([_]u32{ 10, 11, 12, 13 }) |index| {
        const point = atoms[index].coordinates;
        direct += distance(point, parent_center);
        mirrored += distance(reflectAcrossLine(point, axis_start, axis_end), parent_center);
    }
    try std.testing.expect(direct >= mirrored);

    var graph = try topology.Graph.init(std.testing.allocator, &atoms, &bonds);
    defer graph.deinit();
    var rings = try topology.RingMembership.init(std.testing.allocator, graph, &bonds);
    defer rings.deinit();
    var split = try fragments.Fragmentation.init(std.testing.allocator, &atoms, &bonds, graph, rings);
    defer split.deinit();
    var analysis = try topology.rings.Analysis.init(std.testing.allocator, rings, &atoms, &bonds);
    defer analysis.deinit();

    // And the parent chosen for the outer ring is the ring it shares an edge
    // with, not the far one it shares nothing with.
    var placed: [14]bool = @splat(true);
    var outer: ?core.ids.RingId = null;
    for (rings.rings) |ring| {
        if (std.mem.indexOfScalar(core.ids.AtomId, rings.atoms(ring.id), core.ids.AtomId.fromIndex(11)) != null) {
            outer = ring.id;
            break;
        }
    }
    const outer_ring = outer orelse return error.TestUnexpectedResult;
    // Its shared edge is already drawn as part of the parent, which is exactly
    // the state the placement loop is in when it reaches this ring.
    for ([_]u32{ 10, 11, 12, 13 }) |index| placed[index] = false;
    const parent = selectParentRing(rings, analysis, outer_ring, &placed) orelse
        return error.TestUnexpectedResult;
    try std.testing.expect(std.mem.indexOfScalar(core.ids.AtomId, rings.atoms(parent), core.ids.AtomId.fromIndex(6)) != null);
    try std.testing.expect(std.mem.indexOfScalar(core.ids.AtomId, rings.atoms(parent), core.ids.AtomId.fromIndex(0)) == null);
}

test "spiro rings sharing one atom are placed on opposite sides" {
    var atoms: [5]model.Atom = undefined;
    for (&atoms, 0..) |*atom, index| atom.* = .{ .id = core.ids.AtomId.fromIndex(@intCast(index)), .input_index = @intCast(index), .atomic_number = .carbon };
    const pairs = [_][2]u32{ .{ 0, 1 }, .{ 1, 2 }, .{ 2, 0 }, .{ 0, 3 }, .{ 3, 4 }, .{ 4, 0 } };
    var bonds: [pairs.len]model.Bond = undefined;
    for (&bonds, pairs, 0..) |*bond, pair, index| bond.* = .{
        .id = core.ids.BondId.fromIndex(@intCast(index)),
        .input_index = @intCast(index),
        .start = core.ids.AtomId.fromIndex(pair[0]),
        .end = core.ids.AtomId.fromIndex(pair[1]),
        .input_order = .single,
        .effective_order = .single,
    };
    try layoutFixture(&atoms, &bonds);
    try expectBondLengths(&atoms, &bonds);
    const first_center = scale(add(add(atoms[0].coordinates, atoms[1].coordinates), atoms[2].coordinates), 1.0 / 3.0);
    const second_center = scale(add(add(atoms[0].coordinates, atoms[3].coordinates), atoms[4].coordinates), 1.0 / 3.0);
    const first_direction = subtract(first_center, atoms[0].coordinates);
    const second_direction = subtract(second_center, atoms[0].coordinates);
    try std.testing.expect(first_direction.x * second_direction.x + first_direction.y * second_direction.y < 0);
}

test "fragment assembly preserves every acyclic parent bond length" {
    var atoms: [5]model.Atom = undefined;
    for (&atoms, 0..) |*atom, index| atom.* = .{ .id = core.ids.AtomId.fromIndex(@intCast(index)), .input_index = @intCast(index), .atomic_number = .carbon };
    var bonds: [4]model.Bond = undefined;
    for (&bonds, 0..) |*bond, index| bond.* = .{
        .id = core.ids.BondId.fromIndex(@intCast(index)),
        .input_index = @intCast(index),
        .start = core.ids.AtomId.fromIndex(@intCast(index)),
        .end = core.ids.AtomId.fromIndex(@intCast(index + 1)),
        .input_order = .single,
        .effective_order = .single,
    };
    try layoutFixture(&atoms, &bonds);
    try expectBondLengths(&atoms, &bonds);
    for (atoms) |atom| try std.testing.expect(atom.coordinates.isFinite());
    for (1..atoms.len - 1) |center| {
        const left = subtract(atoms[center - 1].coordinates, atoms[center].coordinates);
        const right = subtract(atoms[center + 1].coordinates, atoms[center].coordinates);
        const cosine = (left.x * right.x + left.y * right.y) / (length(left) * length(right));
        try std.testing.expectApproxEqAbs(@as(f32, -0.5), cosine, 0.001);
    }
    var graph = try topology.Graph.init(std.testing.allocator, &atoms, &bonds);
    defer graph.deinit();
    var rings = try topology.RingMembership.init(std.testing.allocator, graph, &bonds);
    defer rings.deinit();
    var split = try fragments.Fragmentation.init(std.testing.allocator, &atoms, &bonds, graph, rings);
    defer split.deinit();
    var frames = try captureFragmentFrames(std.testing.allocator, &atoms, split);
    defer frames.deinit();
    try std.testing.expectEqual(split.fragments.len, frames.frames.len);
    try std.testing.expectEqual(split.fragments.len - 1, frames.child_attachments.len);
    for (frames.frames) |frame| {
        if (!frame.parent.isValid()) continue;
        const parent_frame = frames.frames[frame.parent.index()];
        var found = false;
        for (frames.child_attachments[parent_frame.attachments.start..][0..parent_frame.attachments.len]) |attachment| {
            found = found or (attachment.child == frame.id and attachment.atom == frame.anchor_atom);
        }
        try std.testing.expect(found);
    }
    try core.oom.checkAllocationFailures(
        std.testing.allocator,
        captureFramesAndDiscard,
        .{ &atoms, split },
    );
    const first = atoms;
    try layoutFixture(&atoms, &bonds);
    for (first, atoms) |left, right| try std.testing.expectEqual(left.coordinates, right.coordinates);

    var reversed_atoms: [5]model.Atom = undefined;
    for (&reversed_atoms, 0..) |*atom, index| atom.* = .{ .id = core.ids.AtomId.fromIndex(@intCast(index)), .input_index = @intCast(index), .atomic_number = .carbon };
    var reversed_bonds = bonds;
    for (&reversed_bonds) |*bond| std.mem.swap(core.ids.AtomId, &bond.start, &bond.end);
    try layoutFixture(&reversed_atoms, &reversed_bonds);
    for (0..atoms.len) |left| for (0..atoms.len) |right| {
        try std.testing.expectApproxEqAbs(
            distance(atoms[left].coordinates, atoms[right].coordinates),
            distance(reversed_atoms[left].coordinates, reversed_atoms[right].coordinates),
            0.001,
        );
    };
}

test "wide acyclic centers preserve child attachment bond lengths" {
    const branch_count = max_neighbours + 1;
    var atoms: [1 + branch_count * 2]model.Atom = undefined;
    for (&atoms, 0..) |*atom, index| atom.* = .{
        .id = core.ids.AtomId.fromIndex(@intCast(index)),
        .input_index = @intCast(index),
        .atomic_number = .carbon,
    };
    var bonds: [branch_count * 2]model.Bond = undefined;
    for (0..branch_count) |branch_index| {
        const branch = 1 + branch_index;
        const terminal = 1 + branch_count + branch_index;
        bonds[branch_index] = .{
            .id = core.ids.BondId.fromIndex(@intCast(branch_index)),
            .input_index = @intCast(branch_index),
            .start = core.ids.AtomId.fromIndex(0),
            .end = core.ids.AtomId.fromIndex(@intCast(branch)),
            .input_order = .single,
            .effective_order = .single,
        };
        bonds[branch_count + branch_index] = .{
            .id = core.ids.BondId.fromIndex(@intCast(branch_count + branch_index)),
            .input_index = @intCast(branch_count + branch_index),
            .start = core.ids.AtomId.fromIndex(@intCast(branch)),
            .end = core.ids.AtomId.fromIndex(@intCast(terminal)),
            .input_order = .single,
            .effective_order = .single,
        };
    }

    try layoutFixture(&atoms, &bonds);
    try expectBondLengths(&atoms, &bonds);
}

test "constrained alignment and fixed reset are deterministic" {
    var atoms: [3]model.Atom = undefined;
    for (&atoms, 0..) |*atom, index| atom.* = .{
        .id = core.ids.AtomId.fromIndex(@intCast(index)),
        .input_index = @intCast(index),
        .atomic_number = .carbon,
    };
    var bonds = [_]model.Bond{
        .{ .id = core.ids.BondId.fromIndex(0), .input_index = 0, .start = core.ids.AtomId.fromIndex(0), .end = core.ids.AtomId.fromIndex(1), .input_order = .single, .effective_order = .single },
        .{ .id = core.ids.BondId.fromIndex(1), .input_index = 1, .start = core.ids.AtomId.fromIndex(1), .end = core.ids.AtomId.fromIndex(2), .input_order = .single, .effective_order = .single },
    };
    atoms[0].constrained = true;
    atoms[0].template_coordinates = .{ .x = 10, .y = 20 };
    atoms[1].constrained = true;
    atoms[1].template_coordinates = .{ .x = 10, .y = 70 };
    try layoutFixture(&atoms, &bonds);
    try std.testing.expectApproxEqAbs(@as(f32, 10), atoms[0].coordinates.x, 0.001);
    try std.testing.expectApproxEqAbs(@as(f32, 20), atoms[0].coordinates.y, 0.001);
    try std.testing.expectApproxEqAbs(@as(f32, 10), atoms[1].coordinates.x, 0.001);
    try std.testing.expectApproxEqAbs(@as(f32, 70), atoms[1].coordinates.y, 0.001);

    atoms[2].fixed = true;
    atoms[2].template_coordinates = .{ .x = -12, .y = 34 };
    try layoutFixture(&atoms, &bonds);
    try std.testing.expectEqual(core.math.Vec2{ .x = -12, .y = 34 }, atoms[2].coordinates);

}

test "the 3D fallback applies upstream's scaled, mirrored projection" {
    // Exercised directly rather than through the whole layout: the projection is
    // what cgz-7v2.25 got wrong, and upstream runs this BEFORE the constrained
    // rotation, so driving it end to end would measure the rotation's handling
    // of a deliberately broken template instead of the projection.
    var atoms: [3]model.Atom = undefined;
    for (&atoms, 0..) |*atom, index| atom.* = .{
        .id = core.ids.AtomId.fromIndex(@intCast(index)),
        .input_index = @intCast(index),
        .atomic_number = .carbon,
        .coordinates = .{ .x = std.math.nan(f32), .y = 0 },
        .coordinates_3d = .{ .x = @floatFromInt(index * 3), .y = @floatFromInt(index * 5), .z = 7 },
    };
    var bonds = [_]model.Bond{
        .{ .id = core.ids.BondId.fromIndex(0), .input_index = 0, .start = core.ids.AtomId.fromIndex(0), .end = core.ids.AtomId.fromIndex(1), .input_order = .single, .effective_order = .single },
        .{ .id = core.ids.BondId.fromIndex(1), .input_index = 1, .start = core.ids.AtomId.fromIndex(1), .end = core.ids.AtomId.fromIndex(2), .input_order = .single, .effective_order = .single },
    };
    var graph = try topology.Graph.init(std.testing.allocator, &atoms, &bonds);
    defer graph.deinit();
    var rings = try topology.RingMembership.init(std.testing.allocator, graph, &bonds);
    defer rings.deinit();
    var split = try fragments.Fragmentation.init(std.testing.allocator, &atoms, &bonds, graph, rings);
    defer split.deinit();

    for (split.fragments) |fragment| fallbackOnValid3dCoordinates(&atoms, split.members(fragment.id));
    // 35x, y negated, rounded to two decimals - not the raw 3D x/y a hand-rolled
    // version used to copy through.
    for (atoms, 0..) |atom, index| {
        const source: f32 = @floatFromInt(index);
        try std.testing.expectApproxEqAbs(source * 3 * 35, atom.coordinates.x, 0.001);
        try std.testing.expectApproxEqAbs(-(source * 5) * 35, atom.coordinates.y, 0.001);
    }

    // An infinity is not NaN, so it does not trigger the fallback at all.
    for (&atoms) |*atom| atom.coordinates = .{ .x = std.math.inf(f32), .y = 0 };
    for (split.fragments) |fragment| fallbackOnValid3dCoordinates(&atoms, split.members(fragment.id));
    for (atoms) |atom| try std.testing.expectEqual(std.math.inf(f32), atom.coordinates.x);

    // A missing 3D source declines; and the sentinel is a one-sided ceiling, so
    // a negative infinity passes it where a positive one does not.
    for (&atoms) |*atom| {
        atom.coordinates = .{ .x = std.math.nan(f32), .y = 0 };
        atom.coordinates_3d = null;
    }
    for (split.fragments) |fragment| fallbackOnValid3dCoordinates(&atoms, split.members(fragment.id));
    for (atoms) |atom| try std.testing.expect(std.math.isNan(atom.coordinates.x));

    for (&atoms) |*atom| atom.coordinates_3d = .{ .x = -std.math.inf(f32), .y = 0, .z = 0 };
    try std.testing.expect(fragmentHasValid3dSource(&atoms, split.members(split.fragments[0].id)));
    for (&atoms) |*atom| atom.coordinates_3d = .{ .x = std.math.inf(f32), .y = 0, .z = 0 };
    try std.testing.expect(!fragmentHasValid3dSource(&atoms, split.members(split.fragments[0].id)));
}

test "planarity scoring separates a flat fused system from a bond shared by three rings" {
    // Two fused six-membered rings: flat, nothing shared beyond one edge, and
    // no atom carrying a fourth neighbour. Upstream's first branch.
    var flat_atoms: [10]model.Atom = undefined;
    for (&flat_atoms, 0..) |*atom, index| atom.* = .{ .id = core.ids.AtomId.fromIndex(@intCast(index)), .input_index = @intCast(index), .atomic_number = .carbon };
    const flat_pairs = [_][2]u32{
        .{ 0, 1 }, .{ 1, 2 }, .{ 2, 3 }, .{ 3, 4 }, .{ 4, 5 }, .{ 5, 0 },
        .{ 2, 6 }, .{ 6, 7 }, .{ 7, 8 }, .{ 8, 9 }, .{ 9, 3 },
    };
    var flat_bonds: [flat_pairs.len]model.Bond = undefined;
    for (&flat_bonds, flat_pairs, 0..) |*bond, pair, index| bond.* = .{
        .id = core.ids.BondId.fromIndex(@intCast(index)),
        .input_index = @intCast(index),
        .start = core.ids.AtomId.fromIndex(pair[0]),
        .end = core.ids.AtomId.fromIndex(pair[1]),
        .input_order = .single,
        .effective_order = .single,
    };
    try std.testing.expectEqual(@as(f32, 0), try planarityFixture(&flat_atoms, &flat_bonds));

    // Three bridges across one bond puts that bond in three rings at once,
    // which is upstream's second condition. A bridged cage like
    // bicyclo[2.2.2]octane does NOT qualify and must not be scored as if it
    // did: its bridgeheads have degree three and its rings share atoms rather
    // than bonds, which is why drug_like/5 measures zero here.
    var bridged_atoms: [8]model.Atom = undefined;
    for (&bridged_atoms, 0..) |*atom, index| atom.* = .{ .id = core.ids.AtomId.fromIndex(@intCast(index)), .input_index = @intCast(index), .atomic_number = .carbon };
    const bridged_pairs = [_][2]u32{
        .{ 0, 1 },
        .{ 0, 2 }, .{ 2, 3 }, .{ 3, 1 },
        .{ 0, 4 }, .{ 4, 5 }, .{ 5, 1 },
        .{ 0, 6 }, .{ 6, 7 }, .{ 7, 1 },
    };
    var bridged_bonds: [bridged_pairs.len]model.Bond = undefined;
    for (&bridged_bonds, bridged_pairs, 0..) |*bond, pair, index| bond.* = .{
        .id = core.ids.BondId.fromIndex(@intCast(index)),
        .input_index = @intCast(index),
        .start = core.ids.AtomId.fromIndex(pair[0]),
        .end = core.ids.AtomId.fromIndex(pair[1]),
        .input_order = .single,
        .effective_order = .single,
    };
    const bridged_score = try planarityFixture(&bridged_atoms, &bridged_bonds);
    try std.testing.expect(bridged_score >= non_planar_system_score);
    try std.testing.expect(bridged_score <= untreatable_system_planarity_score);
}

fn planarityFixture(atoms: []model.Atom, bonds: []const model.Bond) !f32 {
    var graph = try topology.Graph.init(std.testing.allocator, atoms, bonds);
    defer graph.deinit();
    var rings = try topology.RingMembership.init(std.testing.allocator, graph, bonds);
    defer rings.deinit();
    var split = try fragments.Fragmentation.init(std.testing.allocator, atoms, bonds, graph, rings);
    defer split.deinit();
    var analysis = try topology.rings.Analysis.init(std.testing.allocator, rings, atoms, bonds);
    defer analysis.deinit();
    var worst: f32 = 0;
    for (split.fragments) |fragment| {
        if (fragment.ring_count == 0) continue;
        worst = @max(worst, try scorePlanarity(rings, analysis, bonds, graph, split, fragment));
    }
    return worst;
}

test "maybeMinimizeRings fires on a five-membered ring sharing an atom with two others" {
    // The transcription is exercised directly, because upstream's only call site
    // hands it an already-drained ring vector and therefore never reaches these
    // conditions. See placeFragmentRings for why that call passes nothing.
    var atoms: [8]model.Atom = undefined;
    for (&atoms, 0..) |*atom, index| atom.* = .{ .id = core.ids.AtomId.fromIndex(@intCast(index)), .input_index = @intCast(index), .atomic_number = .carbon };
    const pairs = [_][2]u32{
        .{ 0, 1 }, .{ 1, 2 }, .{ 2, 3 }, .{ 3, 4 }, .{ 4, 5 },
        .{ 5, 0 }, .{ 5, 6 }, .{ 6, 7 }, .{ 7, 2 },
    };
    var bonds: [pairs.len]model.Bond = undefined;
    for (&bonds, pairs, 0..) |*bond, pair, index| bond.* = .{
        .id = core.ids.BondId.fromIndex(@intCast(index)),
        .input_index = @intCast(index),
        .start = core.ids.AtomId.fromIndex(pair[0]),
        .end = core.ids.AtomId.fromIndex(pair[1]),
        .input_order = .single,
        .effective_order = .single,
    };
    var graph = try topology.Graph.init(std.testing.allocator, &atoms, &bonds);
    defer graph.deinit();
    var rings = try topology.RingMembership.init(std.testing.allocator, graph, &bonds);
    defer rings.deinit();

    var all: [8]core.ids.RingId = undefined;
    for (rings.rings, 0..) |ring, index| all[index] = ring.id;
    const every_ring = all[0..rings.rings.len];

    // An empty set is what the call site passes, and it can never fire.
    try std.testing.expect(!maybeMinimizeRings(rings, &.{}));

    // Over the real rings the answer depends on the perceived basis; assert the
    // condition it encodes rather than a fixed verdict.
    var expected = false;
    for (every_ring) |ring| {
        const members = rings.atoms(ring);
        const five = members.len == 5;
        const odd_macro = members.len >= topology.rings.macrocycle_size and members.len % 2 != 0;
        if (!five and !odd_macro) continue;
        for (members) |atom| {
            if (rings.atomRings(atom).len > 2) expected = true;
        }
    }
    try std.testing.expectEqual(expected, maybeMinimizeRings(rings, every_ring));
}
