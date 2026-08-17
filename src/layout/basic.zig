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
            const start = members[0];
            atoms[start.index()].coordinates = .{};
            placed[start.index()] = true;
            queue[0] = start;
            tail = 1;
        }
        while (head < tail) : (head += 1) {
            const center = queue[head];
            var unplaced_count: usize = 0;
            for (graph.neighbors(center)) |neighbor| {
                if (fragmentation.atom_fragment[neighbor.index()] == fragment.id and !placed[neighbor.index()]) unplaced_count += 1;
            }
            if (unplaced_count == 0) continue;
            const base_angle = parentAngle(atoms, graph, fragmentation, fragment.id, center, placed);
            var generated: usize = 0;
            for (graph.neighbors(center)) |neighbor| {
                if (fragmentation.atom_fragment[neighbor.index()] != fragment.id or placed[neighbor.index()]) continue;
                const spread = if (unplaced_count == 1)
                    @as(f32, 0)
                else
                    (@as(f32, @floatFromInt(generated)) - @as(f32, @floatFromInt(unplaced_count - 1)) * 0.5) * (2 * pi / 3);
                const angle = base_angle + spread;
                atoms[neighbor.index()].coordinates = .{
                    .x = atoms[center.index()].coordinates.x + @cos(angle) * bond_length,
                    .y = atoms[center.index()].coordinates.y + @sin(angle) * bond_length,
                };
                placed[neighbor.index()] = true;
                queue[tail] = neighbor;
                tail += 1;
                generated += 1;
            }
        }
    }
    try assembleFragments(atoms, bonds, graph, fragmentation);
    alignConstrainedMainFragments(atoms, fragmentation);
    fallbackOn3dCoordinates(atoms, fragmentation);
    restoreFixedCoordinates(atoms);
    return outcome;
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
            if (try macrocycle.generateRingShape(
                allocator,
                ring.id,
                atoms,
                bonds,
                graph,
                membership,
                placed,
                options.force_open_macrocycles,
            ) == .matched) {
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
            try alignFusedRing(atoms, ordered, local, placed);
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
    return maybeMinimizeRings(membership, fragmentation, fragment);
}

/// Upstream's CoordgenMinimizer::maybeMinimizeRings, evaluated over the rings
/// of one fragment: a five-membered ring, or an odd-sized macrocycle, holding
/// an atom that belongs to more than two rings requires minimization. It is
/// called only for systems placed without a template, which is why it lives at
/// the end of the fallback path rather than at the top.
fn maybeMinimizeRings(
    membership: topology.RingMembership,
    fragmentation: fragments.Fragmentation,
    fragment: fragments.Fragment,
) bool {
    for (membership.rings) |ring| {
        const ring_atoms = membership.atoms(ring.id);
        if (ring_atoms.len == 0) continue;
        if (fragmentation.atom_fragment[ring_atoms[0].index()] != fragment.id) continue;
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

fn alignFusedRing(atoms: []model.Atom, ordered: []const core.ids.AtomId, local: []const core.math.Vec2, placed: []bool) core.errors.Error!void {
    var first: ?usize = null;
    var last: usize = 0;
    for (ordered, 0..) |atom, index| if (placed[atom.index()]) {
        if (first == null) first = index;
        last = index;
    };
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
    const existing_center = placedCenter(atoms, placed);
    for (ordered, local) |atom, coordinate| {
        if (placed[atom.index()]) continue;
        const candidate = transformFromPivot(coordinate, source_first, target_first, rotation);
        const mirror = reflectAcrossLine(candidate, target_first, target_last);
        first_score += distance(candidate, existing_center);
        mirror_score += distance(mirror, existing_center);
    }
    const use_mirror = mirror_score > first_score;
    for (ordered, local) |atom, coordinate| if (!placed[atom.index()]) {
        const candidate = transformFromPivot(coordinate, source_first, target_first, rotation);
        atoms[atom.index()].coordinates = if (use_mirror) reflectAcrossLine(candidate, target_first, target_last) else candidate;
        placed[atom.index()] = true;
    };
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
    bonds: []const model.Bond,
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
            if (!fragment.bond_to_parent.isValid() or fragment.bond_to_parent.index() >= bonds.len) return error.InvalidMapping;
            const bond = bonds[fragment.bond_to_parent.index()];
            if (fragmentation.atom_fragment[bond.start.index()] == fragment.id) {
                frame.anchor_atom = bond.start;
                frame.parent_atom = bond.end;
            } else if (fragmentation.atom_fragment[bond.end.index()] == fragment.id) {
                frame.anchor_atom = bond.end;
                frame.parent_atom = bond.start;
            } else return error.InvalidMapping;
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

fn assembleFragments(atoms: []model.Atom, bonds: []const model.Bond, graph: topology.Graph, fragmentation: fragments.Fragmentation) core.errors.Error!void {
    for (1..fragmentation.fragments.len) |depth| {
        for (fragmentation.fragments) |fragment| {
            if (!fragment.parent.isValid() or fragmentDepth(fragmentation.fragments, fragment.id) != depth) continue;
            if (!fragment.bond_to_parent.isValid()) return error.InvalidMapping;
            const bond = bonds[fragment.bond_to_parent.index()];
            const child_endpoint = if (fragmentation.atom_fragment[bond.start.index()] == fragment.id) bond.start else bond.end;
            const parent_endpoint = if (child_endpoint == bond.start) bond.end else bond.start;
            const child_center = fragmentCenter(atoms, fragmentation.members(fragment.id));
            var outward = openValenceDirection(atoms, graph, fragmentation, fragment.id, parent_endpoint, depth);
            outward = scale(outward, bond_length / length(outward));
            const target = add(atoms[parent_endpoint.index()].coordinates, outward);
            var child_direction = subtract(child_center, atoms[child_endpoint.index()].coordinates);
            if (length(child_direction) < 0.0001) child_direction = .{ .x = 1 };
            var desired_child_direction = outward;
            if (graph.degree(child_endpoint) == 2) {
                desired_child_direction = rotateVector(scale(outward, -1), 2 * pi / 3);
            }
            const angle = std.math.atan2(desired_child_direction.y, desired_child_direction.x) - std.math.atan2(child_direction.y, child_direction.x);
            const source = atoms[child_endpoint.index()].coordinates;
            for (fragmentation.members(fragment.id)) |atom| {
                atoms[atom.index()].coordinates = transformFromPivot(atoms[atom.index()].coordinates, source, target, angle);
            }
        }
    }
}

fn openValenceDirection(
    atoms: []const model.Atom,
    graph: topology.Graph,
    fragmentation: fragments.Fragmentation,
    child: core.ids.FragmentId,
    center: core.ids.AtomId,
    child_depth: usize,
) core.math.Vec2 {
    var direction_sum: core.math.Vec2 = .{};
    var first_direction: core.math.Vec2 = .{};
    var count: usize = 0;
    for (graph.neighbors(center)) |neighbor| {
        const neighbor_fragment = fragmentation.atom_fragment[neighbor.index()];
        if (neighbor_fragment == child or !neighbor_fragment.isValid() or
            fragmentDepth(fragmentation.fragments, neighbor_fragment) >= child_depth) continue;
        const direction = subtract(atoms[neighbor.index()].coordinates, atoms[center.index()].coordinates);
        if (length(direction) < 0.0001) continue;
        const normalized = scale(direction, 1 / length(direction));
        if (count == 0) first_direction = normalized;
        direction_sum = add(direction_sum, normalized);
        count += 1;
    }
    if (count == 0) return .{ .x = 1 };
    if (count == 1) return rotateVector(first_direction, -2 * pi / 3);
    if (length(direction_sum) < 0.0001) return rotateVector(first_direction, pi);
    return scale(direction_sum, -1 / length(direction_sum));
}

fn rotateVector(value: core.math.Vec2, angle: f32) core.math.Vec2 {
    const cosine = @cos(angle);
    const sine = @sin(angle);
    return .{ .x = value.x * cosine - value.y * sine, .y = value.x * sine + value.y * cosine };
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

fn alignConstrainedMainFragments(atoms: []model.Atom, fragmentation: fragments.Fragmentation) void {
    for (fragmentation.main_fragments) |main| {
        if (!main.isValid()) continue;
        const fragment = fragmentation.fragments[main.index()];
        if (!fragment.flags.constrained or fragment.flags.fixed) continue;
        const members = fragmentation.members(main);
        var source_center: core.math.Vec2 = .{};
        var target_center: core.math.Vec2 = .{};
        var constrained_count: usize = 0;
        for (members) |atom| if (atoms[atom.index()].constrained) {
            const target = atoms[atom.index()].template_coordinates orelse continue;
            source_center = add(source_center, atoms[atom.index()].coordinates);
            target_center = add(target_center, target);
            constrained_count += 1;
        };
        if (constrained_count == 0) continue;
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

fn restoreFixedCoordinates(atoms: []model.Atom) void {
    for (atoms) |*atom| if (atom.fixed) {
        if (atom.template_coordinates) |coordinate| atom.coordinates = coordinate;
    };
}

fn fallbackOn3dCoordinates(atoms: []model.Atom, fragmentation: fragments.Fragmentation) void {
    for (fragmentation.fragments) |fragment| {
        const members = fragmentation.members(fragment.id);
        var generated_finite = true;
        var source_valid = true;
        for (members) |atom| {
            generated_finite = generated_finite and atoms[atom.index()].coordinates.isFinite();
            source_valid = source_valid and atoms[atom.index()].coordinates_3d != null and atoms[atom.index()].coordinates_3d.?.isFinite();
        }
        if (generated_finite or !source_valid) continue;
        for (members) |atom| {
            const source = atoms[atom.index()].coordinates_3d.?;
            atoms[atom.index()].coordinates = .{ .x = source.x, .y = source.y };
        }
    }
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
    try std.testing.checkAllAllocationFailures(
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
    bonds: []const model.Bond,
    split: fragments.Fragmentation,
) !void {
    var frames = try captureFragmentFrames(allocator, atoms, bonds, split);
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
    try std.testing.checkAllAllocationFailures(
        std.testing.allocator,
        layoutAndDiscard,
        .{ &atoms, &bonds, graph, rings, split },
    );
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
    var frames = try captureFragmentFrames(std.testing.allocator, &atoms, &bonds, split);
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
    try std.testing.checkAllAllocationFailures(
        std.testing.allocator,
        captureFramesAndDiscard,
        .{ &atoms, &bonds, split },
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

test "constrained alignment, fixed reset, and valid 3D fallback are deterministic" {
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

    for (&atoms, 0..) |*atom, index| {
        atom.fixed = false;
        atom.constrained = true;
        atom.template_coordinates = .{ .x = std.math.nan(f32), .y = 0 };
        atom.coordinates_3d = .{ .x = @floatFromInt(index * 3), .y = @floatFromInt(index * 5), .z = 7 };
    }
    try layoutFixture(&atoms, &bonds);
    for (atoms, 0..) |atom, index| {
        try std.testing.expectEqual(@as(f32, @floatFromInt(index * 3)), atom.coordinates.x);
        try std.testing.expectEqual(@as(f32, @floatFromInt(index * 5)), atom.coordinates.y);
    }
}
