const std = @import("std");
const core = @import("core");
const model = @import("model");
const geometry = @import("geometry");
const topology = @import("topology");
const layout = @import("layout");
const inversions = @import("inversions.zig");

const Vec2 = core.math.Vec2;
const Bounds = struct { min: Vec2, max: Vec2 };

pub const ProximityRelation = struct {
    start: core.ids.AtomId,
    end: core.ids.AtomId,
};

pub fn collectProximityRelations(
    allocator: std.mem.Allocator,
    working: model.WorkingGraph,
) core.errors.Error![]ProximityRelation {
    var count: usize = 0;
    for (working.proximity_relations) |relation| {
        const endpoints = try proximityEndpoints(
            working.bonds,
            working.extra_bonds,
            working.residue_interactions,
            relation,
        );
        count += @intFromBool(!isResidueRelation(working.residues, relation, endpoints));
    }
    const output = allocator.alloc(ProximityRelation, count) catch return error.OutOfMemory;
    errdefer allocator.free(output);
    var output_index: usize = 0;
    for (working.proximity_relations) |relation| {
        const endpoints = try proximityEndpoints(
            working.bonds,
            working.extra_bonds,
            working.residue_interactions,
            relation,
        );
        if (isResidueRelation(working.residues, relation, endpoints)) continue;
        output[output_index] = endpoints;
        output_index += 1;
    }
    return output;
}

/// Upstream keeps residue interactions out of molecule-proximity placement
/// whenever either endpoint is a residue; the dedicated residue phase owns
/// them. Ordinary zero-order bonds remain proximity relations even if callers
/// also attach residue metadata to an endpoint.
fn isResidueRelation(
    residues: []const model.Residue,
    relation: model.ProximityRelation,
    endpoints: ProximityRelation,
) bool {
    if (relation != .residue_interaction) return false;
    for (residues) |residue| {
        if (residue.atom == endpoints.start or residue.atom == endpoints.end) return true;
    }
    return false;
}

fn proximityEndpoints(
    bonds: []const model.Bond,
    extra_bonds: []const model.ExtraBond,
    residue_interactions: []const model.ResidueInteraction,
    relation: model.ProximityRelation,
) core.errors.Error!ProximityRelation {
    return switch (relation) {
        .bond => |id| if (id.index() < bonds.len)
            .{ .start = bonds[id.index()].start, .end = bonds[id.index()].end }
        else
            error.InvalidMapping,
        .extra_bond => |input_index| blk: {
            for (extra_bonds) |bond| if (bond.input_index == input_index) {
                break :blk .{ .start = bond.start, .end = bond.end };
            };
            return error.InvalidMapping;
        },
        .residue_interaction => |input_index| if (input_index < residue_interactions.len)
            .{ .start = residue_interactions[input_index].start, .end = residue_interactions[input_index].end }
        else
            error.InvalidMapping,
    };
}

/// Upstream chooses the component with the most incident proximity records,
/// then the most atoms, retaining the first component on a complete tie.
pub fn selectProximityCenter(graph: topology.Graph, relations: []const ProximityRelation) core.errors.Error!core.ids.MoleculeId {
    if (graph.component_count == 0 or relations.len == 0) return error.InvalidMapping;
    var best = core.ids.MoleculeId.fromIndex(0);
    var best_relations: usize = 0;
    var best_atoms = graph.componentMembers(best).len;
    for (0..graph.component_count) |raw_index| {
        const component = core.ids.MoleculeId.fromIndex(@intCast(raw_index));
        var relation_count: usize = 0;
        for (relations) |relation| {
            const start_component = graph.component(relation.start) orelse return error.InvalidMapping;
            const end_component = graph.component(relation.end) orelse return error.InvalidMapping;
            relation_count += @intFromBool(start_component == component);
            relation_count += @intFromBool(end_component == component and end_component != start_component);
        }
        const atom_count = graph.componentMembers(component).len;
        if (relation_count > best_relations or
            (relation_count == best_relations and atom_count > best_atoms))
        {
            best = component;
            best_relations = relation_count;
            best_atoms = atom_count;
        }
    }
    return best;
}

const WeightedAngle = struct { weight: f32, angle: f32 };

/// Apply upstream's global bestRotation/maybeFlip rules (sketcherMinimizer.cpp
/// :822 and :612). Both run on every component whose fragments are neither
/// fixed nor constrained, ring-containing ones included: upstream's only ring
/// gate is per-atom and applies to the neighbour-pair angle term alone. The
/// fragment ring and peptide terms need the fragmentation, so a caller that
/// has none — the residue and proximity meta-molecules, which carry no rings —
/// passes null and contributes neither (cgz-r31.1).
pub fn orientComponents(
    allocator: std.mem.Allocator,
    atoms: []model.Atom,
    bonds: []const model.Bond,
    graph: topology.Graph,
    rings: topology.RingMembership,
    fragmentation: ?layout.Fragmentation,
) core.errors.Error!void {
    var analysis: ?topology.rings.Analysis = null;
    defer if (analysis) |*value| value.deinit();
    if (fragmentation != null and rings.rings.len != 0) {
        analysis = try topology.rings.Analysis.init(allocator, rings, atoms, bonds);
    }

    for (0..graph.component_count) |raw_component| {
        const component = core.ids.MoleculeId.fromIndex(@intCast(raw_component));
        const members = graph.componentMembers(component);
        var constrained = false;
        for (members) |atom| {
            constrained = constrained or atoms[atom.index()].fixed or atoms[atom.index()].constrained;
        }
        if (constrained) continue;

        var angles: std.ArrayList(WeightedAngle) = .empty;
        defer angles.deinit(allocator);
        // addBestRotationInfoForPeptides runs before every other term.
        try addPeptideRotationAngles(allocator, &angles, atoms, bonds, graph, members);
        for (members) |atom| {
            // Upstream skips ring atoms here and only here.
            if (rings.atomRings(atom).len != 0) continue;
            if (graph.degree(atom) <= 1) continue;
            const neighbors = graph.neighbors(atom);
            for (neighbors[0 .. neighbors.len - 1], 0..) |first, first_index| {
                for (neighbors[first_index + 1 ..]) |second| {
                    var weight: f32 = 6;
                    if (graph.degree(first) != 1) weight += 2;
                    if (graph.degree(second) != 1) weight += 2;
                    // Preserve pinned upstream's duplicated second-neighbor
                    // carbon bonus.
                    if (atoms[second.index()].atomic_number == .carbon) weight += 1;
                    if (atoms[second.index()].atomic_number == .carbon) weight += 1;
                    if (atoms[first.index()].formal_charge == 0) weight += 1;
                    if (atoms[second.index()].formal_charge == 0) weight += 1;
                    const direction = geometry.subtract(atoms[first.index()].coordinates, atoms[second.index()].coordinates);
                    try addAngle(allocator, &angles, weight, std.math.atan2(-direction.y, direction.x));
                }
            }
        }
        for (bonds) |bond| {
            if (graph.component(bond.start) != component or bond.effective_order == .zero or bond.skip) continue;
            const direction = geometry.subtract(atoms[bond.end.index()].coordinates, atoms[bond.start.index()].coordinates);
            var angle = roundToTwo(std.math.atan2(-direction.y, direction.x));
            while (angle <= 0) angle += std.math.pi;
            for (0..6) |index| {
                var weight: f32 = if (index == 1 or index == 5) 5 else if (index == 0 or index == 3) 1.5 else 1;
                if (bond.effective_order == .double and index == 3 and
                    (graph.degree(bond.start) == 1 or graph.degree(bond.end) == 1)) weight += 1.5;
                if (graph.degree(bond.start) == 1 and graph.degree(bond.end) == 1 and index == 0) weight += 10;
                try addAngle(allocator, &angles, weight, angle);
                angle += std.math.pi / 6.0;
                if (angle > std.math.pi) angle -= std.math.pi;
            }
        }
        if (fragmentation) |parts| if (analysis) |fusion| {
            try addFragmentRingAngles(allocator, &angles, atoms, graph, rings, fusion, parts, component);
        };
        if (angles.items.len > 1 and
            angles.items[angles.items.len - 1].angle - angles.items[0].angle >= std.math.pi - 2 * geometry.epsilon)
        {
            angles.items[0].weight += angles.items[angles.items.len - 1].weight;
            _ = angles.pop();
        }
        if (angles.items.len != 0) {
            var best: usize = 0;
            for (angles.items, 0..) |candidate, index| if (candidate.weight > angles.items[best].weight) {
                best = index;
            };
            const center = componentCenter(atoms, members);
            const sine = -std.math.sin(angles.items[best].angle);
            const cosine = std.math.cos(angles.items[best].angle);
            for (members) |atom| {
                const relative = geometry.subtract(atoms[atom.index()].coordinates, center);
                atoms[atom.index()].coordinates = geometry.add(center, geometry.rotate(relative, sine, cosine));
            }
        }
        try maybeFlipComponent(allocator, atoms, bonds, graph, rings, analysis, fragmentation, component, members);
    }
}

/// One fragment's rings, in ring-index order. Upstream assigns a ring to the
/// fragment owning its first atom (CoordgenFragmenter.cpp:72) and appends in
/// molecule ring order, so `rings[0]` and `rings[1]` below are that order.
const FragmentRings = struct {
    ids: []core.ids.RingId,
    len: usize,

    fn slice(self: *const FragmentRings) []const core.ids.RingId {
        return self.ids[0..self.len];
    }
};

/// Walk the component's fragments in id order, handing each one's ring list to
/// `visit`. Threaded rather than materialised so the two callers — bestRotation's
/// angle terms and maybeFlip's score terms — share one grouping.
fn forEachFragmentRings(
    allocator: std.mem.Allocator,
    graph: topology.Graph,
    rings: topology.RingMembership,
    parts: layout.Fragmentation,
    component: core.ids.MoleculeId,
    context: anytype,
    comptime visit: fn (@TypeOf(context), []const core.ids.RingId) core.errors.Error!void,
) core.errors.Error!void {
    if (rings.rings.len == 0) return;
    const owner = allocator.alloc(core.ids.FragmentId, rings.rings.len) catch return error.OutOfMemory;
    defer allocator.free(owner);
    for (owner, 0..) |*value, index| {
        const ring = core.ids.RingId.fromIndex(@intCast(index));
        const members = rings.atoms(ring);
        if (members.len == 0) return error.InvalidMapping;
        value.* = parts.atom_fragment[members[0].index()];
    }
    var buffer: FragmentRings = .{
        .ids = allocator.alloc(core.ids.RingId, rings.rings.len) catch return error.OutOfMemory,
        .len = 0,
    };
    defer allocator.free(buffer.ids);
    for (parts.fragments) |fragment| {
        const fragment_atoms = parts.members(fragment.id);
        if (fragment_atoms.len == 0) continue;
        if (graph.component(fragment_atoms[0]) != component) continue;
        buffer.len = 0;
        for (owner, 0..) |value, index| {
            if (value != fragment.id) continue;
            buffer.ids[buffer.len] = core.ids.RingId.fromIndex(@intCast(index));
            buffer.len += 1;
        }
        try visit(context, buffer.slice());
    }
}

fn ringCenter(atoms: []const model.Atom, members: []const core.ids.AtomId) Vec2 {
    var result: Vec2 = .{};
    for (members) |atom| {
        result.x += atoms[atom.index()].coordinates.x;
        result.y += atoms[atom.index()].coordinates.y;
    }
    result.x /= @floatFromInt(members.len);
    result.y /= @floatFromInt(members.len);
    return result;
}

const RingAngleContext = struct {
    allocator: std.mem.Allocator,
    angles: *std.ArrayList(WeightedAngle),
    atoms: []const model.Atom,
    rings: topology.RingMembership,
    fusion: topology.rings.Analysis,
};

/// bestRotation's per-fragment ring terms (sketcherMinimizer.cpp:910-975).
fn addFragmentRingAngles(
    allocator: std.mem.Allocator,
    angles: *std.ArrayList(WeightedAngle),
    atoms: []const model.Atom,
    graph: topology.Graph,
    rings: topology.RingMembership,
    fusion: topology.rings.Analysis,
    parts: layout.Fragmentation,
    component: core.ids.MoleculeId,
) core.errors.Error!void {
    var context = RingAngleContext{
        .allocator = allocator,
        .angles = angles,
        .atoms = atoms,
        .rings = rings,
        .fusion = fusion,
    };
    try forEachFragmentRings(allocator, graph, rings, parts, component, &context, addRingAnglesOfFragment);
}

fn addRingAnglesOfFragment(
    context: *RingAngleContext,
    fragment_rings: []const core.ids.RingId,
) core.errors.Error!void {
    if (fragment_rings.len == 2) {
        const first = ringCenter(context.atoms, context.rings.atoms(fragment_rings[0]));
        const second = ringCenter(context.atoms, context.rings.atoms(fragment_rings[1]));
        const direction = geometry.subtract(second, first);
        try addAngle(context.allocator, context.angles, 25, std.math.atan2(-direction.y, direction.x));
        return;
    }
    if (fragment_rings.len == 3) {
        var first = ringCenter(context.atoms, context.rings.atoms(fragment_rings[0]));
        var second = ringCenter(context.atoms, context.rings.atoms(fragment_rings[1]));
        // The first ring fused to exactly two others across two two-atom
        // fusions replaces both reference points with the fusion midpoints.
        for (fragment_rings) |ring| {
            const fusions = context.fusion.fusedWith(ring);
            if (fusions.len != 2) continue;
            const left = context.fusion.fusionAtoms(fusions[0]);
            const right = context.fusion.fusionAtoms(fusions[1]);
            if (left.len != 2 or right.len != 2) continue;
            first = fusionMidpoint(context.atoms, left[0], left[1]);
            second = fusionMidpoint(context.atoms, right[0], right[1]);
            break;
        }
        const direction = geometry.subtract(second, first);
        try addAngle(context.allocator, context.angles, 50, std.math.atan2(-direction.y, direction.x));
        return;
    }
    for (fragment_rings) |ring| {
        if (context.rings.atoms(ring).len != 6) continue;
        for (context.fusion.fusedWith(ring)) |entry| {
            const shared = context.fusion.fusionAtoms(entry);
            if (shared.len != 2) continue;
            const direction = geometry.subtract(
                context.atoms[shared[0].index()].coordinates,
                context.atoms[shared[1].index()].coordinates,
            );
            try addAngle(
                context.allocator,
                context.angles,
                25,
                std.math.atan2(-direction.y, direction.x) - std.math.pi * 0.5,
            );
        }
    }
}

fn fusionMidpoint(atoms: []const model.Atom, left: core.ids.AtomId, right: core.ids.AtomId) Vec2 {
    return .{
        .x = (atoms[left.index()].coordinates.x + atoms[right.index()].coordinates.x) * 0.5,
        .y = (atoms[left.index()].coordinates.y + atoms[right.index()].coordinates.y) * 0.5,
    };
}

/// The direction each alpha carbon's amino nitrogen takes from its cheto
/// carbon. Upstream reads it twice — once as a rotation angle at weight 1000
/// and once as a flip score of +/-100 — from the same classification, whose
/// two-of-each-class guard lives in the peptide *constraint* builder and not in
/// the classifiers, so no count threshold applies here.
fn forEachPeptideDirection(
    atoms: []const model.Atom,
    bonds: []const model.Bond,
    graph: topology.Graph,
    members: []const core.ids.AtomId,
    context: anytype,
    comptime visit: fn (@TypeOf(context), Vec2) core.errors.Error!void,
) core.errors.Error!void {
    for (members) |candidate| {
        if (atoms[candidate.index()].atomic_number != .carbon) continue;
        if (try inversions.isChetoCarbon(candidate, atoms, bonds, graph)) continue;
        var amino: ?core.ids.AtomId = null;
        var cheto: ?core.ids.AtomId = null;
        var bonded_to_cheto = false;
        var bonded_to_amino = false;
        for (graph.neighbors(candidate)) |neighbor| {
            if (atoms[neighbor.index()].atomic_number == .nitrogen) {
                bonded_to_amino = true;
                amino = neighbor;
            } else if (try inversions.isChetoCarbon(neighbor, atoms, bonds, graph)) {
                bonded_to_cheto = true;
                cheto = neighbor;
            }
        }
        if (!bonded_to_cheto or !bonded_to_amino) continue;
        const amino_atom = amino orelse continue;
        const cheto_atom = cheto orelse continue;
        try visit(context, geometry.subtract(
            atoms[amino_atom.index()].coordinates,
            atoms[cheto_atom.index()].coordinates,
        ));
    }
}

const PeptideAngleContext = struct {
    allocator: std.mem.Allocator,
    angles: *std.ArrayList(WeightedAngle),
};

fn addPeptideRotationAngles(
    allocator: std.mem.Allocator,
    angles: *std.ArrayList(WeightedAngle),
    atoms: []const model.Atom,
    bonds: []const model.Bond,
    graph: topology.Graph,
    members: []const core.ids.AtomId,
) core.errors.Error!void {
    var context = PeptideAngleContext{ .allocator = allocator, .angles = angles };
    try forEachPeptideDirection(atoms, bonds, graph, members, &context, addPeptideAngle);
}

fn addPeptideAngle(context: *PeptideAngleContext, direction: Vec2) core.errors.Error!void {
    try addAngle(context.allocator, context.angles, 1000, std.math.atan2(-direction.y, direction.x));
}

fn addPeptideFlipScore(score_x: *f32, direction: Vec2) core.errors.Error!void {
    const peptide_score: f32 = 100;
    if (direction.x > 0) score_x.* -= peptide_score else score_x.* += peptide_score;
}

const RingFlipContext = struct {
    score_x: *f32,
    score_y: *f32,
    atoms: []const model.Atom,
    rings: topology.RingMembership,
    fusion: topology.rings.Analysis,
};

fn addAngle(allocator: std.mem.Allocator, angles: *std.ArrayList(WeightedAngle), weight: f32, raw_angle: f32) core.errors.Error!void {
    var angle = roundToTwo(raw_angle);
    while (angle <= 0) angle += std.math.pi;
    for (angles.items, 0..) |*candidate, index| {
        if (candidate.angle < angle - geometry.epsilon) continue;
        if (candidate.angle - angle < geometry.epsilon and candidate.angle - angle > -geometry.epsilon) {
            candidate.weight += weight;
        } else {
            angles.insert(allocator, index, .{ .weight = weight, .angle = angle }) catch return error.OutOfMemory;
        }
        return;
    }
    angles.append(allocator, .{ .weight = weight, .angle = angle }) catch return error.OutOfMemory;
}

fn maybeFlipComponent(
    allocator: std.mem.Allocator,
    atoms: []model.Atom,
    bonds: []const model.Bond,
    graph: topology.Graph,
    rings: topology.RingMembership,
    analysis: ?topology.rings.Analysis,
    fragmentation: ?layout.Fragmentation,
    component: core.ids.MoleculeId,
    members: []const core.ids.AtomId,
) core.errors.Error!void {
    if (members.len < 2) return;
    var score_x: f32 = 0;
    var score_y: f32 = 0;
    // maybeFlipPeptides runs before every other term and scores X only.
    try forEachPeptideDirection(atoms, bonds, graph, members, &score_x, addPeptideFlipScore);
    if (fragmentation) |parts| if (analysis) |fusion| {
        var context = RingFlipContext{
            .score_x = &score_x,
            .score_y = &score_y,
            .atoms = atoms,
            .rings = rings,
            .fusion = fusion,
        };
        try forEachFragmentRings(allocator, graph, rings, parts, component, &context, scoreRingFlipOfFragment);
    };
    const center = componentCenter(atoms, members);
    const bounds = flipBounds(atoms, members);
    const midpoint = Vec2{ .x = (bounds.max.x + bounds.min.x) * 0.5, .y = (bounds.max.y + bounds.min.y) * 0.5 };
    if (midpoint.x - center.x > geometry.epsilon) score_x -= 0.5 else if (midpoint.x - center.x < -geometry.epsilon) score_x += 0.5;
    if (midpoint.y - center.y > geometry.epsilon) score_y += 0.5 else if (midpoint.y - center.y < -geometry.epsilon) score_y -= 0.5;
    for (bonds) |bond| {
        if (graph.component(bond.start) != component or bond.effective_order != .double or bond.skip) continue;
        var difference: ?f32 = null;
        if (graph.degree(bond.start) == 1 and graph.degree(bond.end) > 1) {
            difference = atoms[bond.start.index()].coordinates.y - atoms[bond.end.index()].coordinates.y;
        } else if (graph.degree(bond.end) == 1 and graph.degree(bond.start) > 1) {
            difference = atoms[bond.end.index()].coordinates.y - atoms[bond.start.index()].coordinates.y;
        }
        if (difference) |value| {
            if (value > geometry.epsilon) score_y += 1 else if (value < -geometry.epsilon) score_y -= 1;
        }
    }
    for (members) |atom| {
        if (score_y < 0) atoms[atom.index()].coordinates.y = -atoms[atom.index()].coordinates.y;
        if (score_x < 0) atoms[atom.index()].coordinates.x = -atoms[atom.index()].coordinates.x;
    }
}

/// maybeFlip's own extent loop, quirk included: upstream chains the maximum
/// test behind `else` (sketcherMinimizer.cpp:693-706), so an atom that sets a
/// new minimum is never considered for the maximum — the first atom can never
/// set one. Seeds are upstream's literals, not the first atom's coordinates.
/// `componentBounds` is a true bounding box and belongs to other phases; this
/// is deliberately not that function.
fn flipBounds(atoms: []const model.Atom, members: []const core.ids.AtomId) Bounds {
    var result = Bounds{
        .min = .{ .x = 9999, .y = 9999 },
        .max = .{ .x = -9999, .y = -9999 },
    };
    for (members) |atom| {
        const coordinate = atoms[atom.index()].coordinates;
        if (coordinate.x < result.min.x) {
            result.min.x = coordinate.x;
        } else if (coordinate.x > result.max.x) {
            result.max.x = coordinate.x;
        }
        if (coordinate.y < result.min.y) {
            result.min.y = coordinate.y;
        } else if (coordinate.y > result.max.y) {
            result.max.y = coordinate.y;
        }
    }
    return result;
}

/// maybeFlip's per-fragment ring terms (sketcherMinimizer.cpp:634-691).
fn scoreRingFlipOfFragment(
    context: *RingFlipContext,
    fragment_rings: []const core.ids.RingId,
) core.errors.Error!void {
    if (fragment_rings.len == 2) {
        var first = fragment_rings[0];
        var second = fragment_rings[1];
        var first_center = ringCenter(context.atoms, context.rings.atoms(first));
        var second_center = ringCenter(context.atoms, context.rings.atoms(second));
        if (first_center.x - second_center.x > geometry.epsilon) {
            std.mem.swap(core.ids.RingId, &first, &second);
            std.mem.swap(Vec2, &first_center, &second_center);
        }
        const first_benzene = context.fusion.flags[first.index()].benzene;
        const second_benzene = context.fusion.flags[second.index()].benzene;
        if (second_benzene and !first_benzene) {
            context.score_x.* -= 20;
        } else if (first_benzene and !second_benzene) {
            context.score_x.* += 20;
        }
    }
    if (fragment_rings.len > 3) {
        var center: Vec2 = .{};
        var weighted: Vec2 = .{};
        var total_atoms: usize = 0;
        var total_rings: usize = 0;
        for (fragment_rings) |ring| {
            const ring_atoms = context.rings.atoms(ring);
            if (ring_atoms.len < 4) continue;
            const ring_center = ringCenter(context.atoms, ring_atoms);
            center = geometry.add(center, ring_center);
            weighted.x += ring_center.x * @as(f32, @floatFromInt(ring_atoms.len));
            weighted.y += ring_center.y * @as(f32, @floatFromInt(ring_atoms.len));
            total_atoms += ring_atoms.len;
            total_rings += 1;
        }
        if (total_rings != 0 and total_atoms != 0) {
            center.x /= @floatFromInt(total_rings);
            center.y /= @floatFromInt(total_rings);
            weighted.x /= @floatFromInt(total_atoms);
            weighted.y /= @floatFromInt(total_atoms);
        }
        if (weighted.y - center.y < -geometry.epsilon) {
            context.score_y.* += 50;
        } else if (weighted.y - center.y > geometry.epsilon) {
            context.score_y.* -= 50;
        }
        if (weighted.x - center.x < -geometry.epsilon) {
            context.score_x.* += 50;
        } else if (weighted.x - center.x > geometry.epsilon) {
            context.score_x.* -= 50;
        }
    }
}

fn roundToTwo(value: f32) f32 {
    return @round(value * 100) / 100;
}

/// Mirror one component across the axis perpendicular to its first crossing
/// pair of already-placeable proximity relations, matching upstream's
/// flipIfCrossingInteractions first-hit behavior.
pub fn flipFirstCrossingRelations(
    atoms: []model.Atom,
    graph: topology.Graph,
    relations: []const ProximityRelation,
    component: core.ids.MoleculeId,
    placed: []const bool,
) core.errors.Error!bool {
    if (placed.len != graph.component_count) return error.InvalidMapping;
    for (relations, 0..) |first, first_index| {
        const first_local = localEndpoint(graph, first, component) orelse continue;
        const first_other = if (first_local == first.start) first.end else first.start;
        const first_other_component = graph.component(first_other) orelse return error.InvalidMapping;
        if (first_other_component == component or !placed[first_other_component.index()]) continue;
        for (relations[first_index + 1 ..]) |second| {
            const second_local = localEndpoint(graph, second, component) orelse continue;
            const second_other = if (second_local == second.start) second.end else second.start;
            const second_other_component = graph.component(second_other) orelse return error.InvalidMapping;
            if (second_other_component == component or !placed[second_other_component.index()]) continue;
            if (geometry.segmentIntersection(
                atoms[first.start.index()].coordinates,
                atoms[first.end.index()].coordinates,
                atoms[second.start.index()].coordinates,
                atoms[second.end.index()].coordinates,
            ) == null) continue;

            const first_point = atoms[first_local.index()].coordinates;
            const second_point = atoms[second_local.index()].coordinates;
            const middle = geometry.scale(geometry.add(first_point, second_point), 0.5);
            const axis = geometry.normalize(geometry.subtract(first_point, second_point));
            for (graph.componentMembers(component)) |atom| {
                const relative = geometry.subtract(atoms[atom.index()].coordinates, middle);
                const parallel = geometry.scale(axis, geometry.dot(axis, relative));
                atoms[atom.index()].coordinates = geometry.subtract(atoms[atom.index()].coordinates, geometry.scale(parallel, 2));
                atoms[atom.index()].coordinates.x = @round(atoms[atom.index()].coordinates.x);
                atoms[atom.index()].coordinates.y = @round(atoms[atom.index()].coordinates.y);
            }
            return true;
        }
    }
    return false;
}

fn localEndpoint(graph: topology.Graph, relation: ProximityRelation, component: core.ids.MoleculeId) ?core.ids.AtomId {
    if (graph.component(relation.start) == component) return relation.start;
    if (graph.component(relation.end) == component) return relation.end;
    return null;
}

/// Place one relation-connected child using upstream's ligand/residue-style
/// interaction-site averages and parent free-valence direction. The caller
/// owns BFS order and marks the central component before invoking this.
pub fn placeProximityChild(
    atoms: []model.Atom,
    graph: topology.Graph,
    rings: topology.RingMembership,
    relations: []const ProximityRelation,
    child: core.ids.MoleculeId,
    parent: core.ids.MoleculeId,
    placed: []bool,
) core.errors.Error!bool {
    if (placed.len != graph.component_count or child == parent or !placed[parent.index()]) return error.InvalidMapping;
    const child_members = graph.componentMembers(child);
    const child_center = componentCenter(atoms, child_members);
    var parent_position: Vec2 = .{};
    var parent_addition: Vec2 = .{};
    var child_position: Vec2 = .{};
    var child_direction: Vec2 = .{};
    var count: u32 = 0;
    for (relations) |relation| {
        const first_component = graph.component(relation.start) orelse return error.InvalidMapping;
        const second_component = graph.component(relation.end) orelse return error.InvalidMapping;
        const parent_atom, const child_atom = if (first_component == parent and second_component == child)
            .{ relation.start, relation.end }
        else if (second_component == parent and first_component == child)
            .{ relation.end, relation.start }
        else
            continue;
        count += 1;
        parent_position = geometry.add(parent_position, atoms[parent_atom.index()].coordinates);
        var addition = singleAdditionVector(atoms, graph, rings, parent_atom);
        addition = geometry.scale(geometry.normalize(addition), core.math.bond_length * 3);
        parent_addition = geometry.add(parent_addition, addition);
        child_position = geometry.add(child_position, atoms[child_atom.index()].coordinates);
        child_direction = geometry.add(child_direction, geometry.subtract(atoms[child_atom.index()].coordinates, child_center));
    }
    if (count == 0) return false;
    const divisor: f32 = @floatFromInt(count);
    parent_position = geometry.scale(parent_position, 1 / divisor);
    parent_addition = geometry.scale(parent_addition, 1 / divisor);
    child_position = geometry.scale(child_position, 1 / divisor);
    child_direction = geometry.scale(child_direction, 1 / divisor);
    var starting_position = geometry.add(parent_position, parent_addition);
    starting_position = findGridPoint(atoms, graph, placed, starting_position, 0, 0, core.math.bond_length * 1.8);
    const desired = geometry.subtract(starting_position, parent_position);
    const opposite_child = geometry.scale(child_direction, -1);
    const angle = geometry.signedAngle(desired, .{}, opposite_child) / 180 * std.math.pi;
    const sine = std.math.sin(angle);
    const cosine = std.math.cos(angle);
    for (child_members) |atom| {
        const relative = geometry.subtract(atoms[atom.index()].coordinates, child_position);
        atoms[atom.index()].coordinates = geometry.add(geometry.rotate(relative, sine, cosine), starting_position);
        atoms[atom.index()].coordinates.x = @round(atoms[atom.index()].coordinates.x);
        atoms[atom.index()].coordinates.y = @round(atoms[atom.index()].coordinates.y);
    }
    placed[child.index()] = true;
    _ = try flipFirstCrossingRelations(atoms, graph, relations, child, placed);
    return true;
}

/// Place an acyclic connected proximity subsystem using upstream's
/// ligand/residue style. Ring-shaped meta-graphs and central components below
/// the pinned eight-atom cutoff belong to the general alignment path.
pub fn arrangeProximityComponents(
    allocator: std.mem.Allocator,
    atoms: []model.Atom,
    graph: topology.Graph,
    rings: topology.RingMembership,
    relations: []const ProximityRelation,
) core.errors.Error!bool {
    if (relations.len == 0) return false;
    const central = try selectProximityCenter(graph, relations);

    const related = allocator.alloc(bool, graph.component_count) catch return error.OutOfMemory;
    defer allocator.free(related);
    @memset(related, false);
    var unique_edges: usize = 0;
    for (relations, 0..) |relation, index| {
        const first = graph.component(relation.start) orelse return error.InvalidMapping;
        const second = graph.component(relation.end) orelse return error.InvalidMapping;
        if (first == second) continue;
        related[first.index()] = true;
        related[second.index()] = true;
        var duplicate = false;
        for (relations[0..index]) |previous| {
            const previous_first = graph.component(previous.start) orelse return error.InvalidMapping;
            const previous_second = graph.component(previous.end) orelse return error.InvalidMapping;
            if ((previous_first == first and previous_second == second) or
                (previous_first == second and previous_second == first))
            {
                duplicate = true;
                break;
            }
        }
        unique_edges += @intFromBool(!duplicate);
    }
    var related_count: usize = 0;
    for (related) |value| related_count += @intFromBool(value);
    if (related_count != graph.component_count) return false;
    if (unique_edges >= related_count or graph.componentMembers(central).len < 8) {
        return try arrangeGeneralProximityComponents(allocator, atoms, graph, rings, relations, unique_edges);
    }

    const placed = allocator.alloc(bool, graph.component_count) catch return error.OutOfMemory;
    defer allocator.free(placed);
    @memset(placed, false);
    placed[central.index()] = true;
    const parents = allocator.alloc(core.ids.MoleculeId, graph.component_count) catch return error.OutOfMemory;
    defer allocator.free(parents);
    @memset(parents, .invalid);
    var queue: std.ArrayList(core.ids.MoleculeId) = .empty;
    defer queue.deinit(allocator);
    queue.append(allocator, central) catch return error.OutOfMemory;
    var head: usize = 0;
    while (head < queue.items.len) : (head += 1) {
        const current = queue.items[head];
        if (current != central) {
            if (placed[current.index()]) continue;
            const parent = parents[current.index()];
            if (!parent.isValid() or !try placeProximityChild(atoms, graph, rings, relations, current, parent, placed)) return error.InvalidMapping;
        }
        for (relations) |relation| {
            const first = graph.component(relation.start) orelse return error.InvalidMapping;
            const second = graph.component(relation.end) orelse return error.InvalidMapping;
            const neighbor = if (first == current and second != current)
                second
            else if (second == current and first != current)
                first
            else
                continue;
            if (placed[neighbor.index()]) continue;
            parents[neighbor.index()] = current;
            queue.append(allocator, neighbor) catch return error.OutOfMemory;
        }
    }
    for (related, placed) |is_related, is_placed| if (is_related and !is_placed) return false;
    return true;
}

/// Run proximity placement on a compact copy of non-residue atoms, preserving
/// full-graph coordinates and IDs outside that active subsystem.
pub fn arrangeProximityComponentsExcluding(
    allocator: std.mem.Allocator,
    atoms: []model.Atom,
    bonds: []const model.Bond,
    relations: []const ProximityRelation,
    excluded_atoms: []const bool,
) core.errors.Error!bool {
    if (excluded_atoms.len != atoms.len) return error.InvalidMapping;
    const old_to_new = allocator.alloc(u32, atoms.len) catch return error.OutOfMemory;
    defer allocator.free(old_to_new);
    @memset(old_to_new, std.math.maxInt(u32));
    var active_count: usize = 0;
    for (excluded_atoms) |excluded| active_count += @intFromBool(!excluded);
    if (active_count == 0) return false;
    const active_atoms = allocator.alloc(model.Atom, active_count) catch return error.OutOfMemory;
    defer allocator.free(active_atoms);
    const new_to_old = allocator.alloc(u32, active_count) catch return error.OutOfMemory;
    defer allocator.free(new_to_old);
    var active_index: usize = 0;
    for (atoms, excluded_atoms, 0..) |atom, excluded, old_index| {
        if (excluded) continue;
        if (active_index >= std.math.maxInt(u32)) return error.TooManyItems;
        old_to_new[old_index] = @intCast(active_index);
        new_to_old[active_index] = @intCast(old_index);
        active_atoms[active_index] = atom;
        active_atoms[active_index].id = core.ids.AtomId.fromIndex(@intCast(active_index));
        active_index += 1;
    }

    var active_bond_count: usize = 0;
    for (bonds) |bond| {
        if (bond.start.index() >= atoms.len or bond.end.index() >= atoms.len) return error.InvalidMapping;
        active_bond_count += @intFromBool(!excluded_atoms[bond.start.index()] and !excluded_atoms[bond.end.index()]);
    }
    const active_bonds = allocator.alloc(model.Bond, active_bond_count) catch return error.OutOfMemory;
    defer allocator.free(active_bonds);
    var bond_index: usize = 0;
    for (bonds) |bond| {
        if (excluded_atoms[bond.start.index()] or excluded_atoms[bond.end.index()]) continue;
        active_bonds[bond_index] = bond;
        active_bonds[bond_index].id = core.ids.BondId.fromIndex(@intCast(bond_index));
        active_bonds[bond_index].start = core.ids.AtomId.fromIndex(old_to_new[bond.start.index()]);
        active_bonds[bond_index].end = core.ids.AtomId.fromIndex(old_to_new[bond.end.index()]);
        bond_index += 1;
    }
    const active_relations = allocator.alloc(ProximityRelation, relations.len) catch return error.OutOfMemory;
    defer allocator.free(active_relations);
    for (relations, active_relations) |relation, *active| {
        if (relation.start.index() >= atoms.len or relation.end.index() >= atoms.len or
            excluded_atoms[relation.start.index()] or excluded_atoms[relation.end.index()]) return error.InvalidMapping;
        active.* = .{
            .start = core.ids.AtomId.fromIndex(old_to_new[relation.start.index()]),
            .end = core.ids.AtomId.fromIndex(old_to_new[relation.end.index()]),
        };
    }

    var graph = try topology.Graph.init(allocator, active_atoms, active_bonds);
    defer graph.deinit();
    var rings = try topology.RingMembership.init(allocator, graph, active_bonds);
    defer rings.deinit();
    if (!try arrangeProximityComponents(allocator, active_atoms, graph, rings, active_relations)) return false;
    for (active_atoms, new_to_old) |atom, old_index| atoms[old_index].coordinates = atom.coordinates;
    return true;
}

/// General meta-graph placement. Upstream recursively lays out one carbon-like
/// meta atom per component, aligns each real component's incident sites to its
/// meta neighbors, then expands that template until 25-unit clashes clear.
fn arrangeGeneralProximityComponents(
    allocator: std.mem.Allocator,
    atoms: []model.Atom,
    graph: topology.Graph,
    rings: topology.RingMembership,
    relations: []const ProximityRelation,
    unique_edge_count: usize,
) core.errors.Error!bool {
    const templates = (try generateMetaCoordinates(allocator, graph, relations, unique_edge_count)) orelse return false;
    defer allocator.free(templates);
    for (0..graph.component_count) |raw_component| {
        const component = core.ids.MoleculeId.fromIndex(@intCast(raw_component));
        if (!try rotateGeneralComponent(allocator, atoms, graph, rings, relations, templates, component)) return false;
    }

    var counter: u32 = 1;
    var clever = false;
    while (true) {
        clever = !clever;
        if (!clever) counter += 1;
        for (0..graph.component_count) |raw_component| {
            const component = core.ids.MoleculeId.fromIndex(@intCast(raw_component));
            const site = try proximitySite(atoms, graph, rings, relations, component);
            const target = geometry.scale(templates[raw_component], @floatFromInt(counter));
            const translation = geometry.subtract(target, site.center);
            for (graph.componentMembers(component)) |atom| {
                atoms[atom.index()].coordinates = geometry.add(atoms[atom.index()].coordinates, translation);
            }
        }
        if (clever) try replaceTerminalSingletons(atoms, graph, rings, relations, counter);
        if (!componentsClash(atoms, graph, core.math.bond_length * 0.5) or counter >= 10) break;
    }
    return true;
}

fn generateMetaCoordinates(
    allocator: std.mem.Allocator,
    graph: topology.Graph,
    relations: []const ProximityRelation,
    unique_edge_count: usize,
) core.errors.Error!?[]Vec2 {
    const meta_atoms = allocator.alloc(model.Atom, graph.component_count) catch return error.OutOfMemory;
    defer allocator.free(meta_atoms);
    for (meta_atoms, 0..) |*atom, index| atom.* = .{
        .id = core.ids.AtomId.fromIndex(@intCast(index)),
        .input_index = @intCast(index),
        .atomic_number = .carbon,
    };
    const meta_bonds = allocator.alloc(model.Bond, unique_edge_count) catch return error.OutOfMemory;
    defer allocator.free(meta_bonds);
    var bond_count: usize = 0;
    for (relations) |relation| {
        const start = graph.component(relation.start) orelse return error.InvalidMapping;
        const end = graph.component(relation.end) orelse return error.InvalidMapping;
        if (start == end) continue;
        var duplicate = false;
        for (meta_bonds[0..bond_count]) |bond| {
            duplicate = duplicate or (bond.start == core.ids.AtomId.fromIndex(start.index()) and bond.end == core.ids.AtomId.fromIndex(end.index())) or
                (bond.start == core.ids.AtomId.fromIndex(end.index()) and bond.end == core.ids.AtomId.fromIndex(start.index()));
        }
        if (duplicate) continue;
        meta_bonds[bond_count] = .{
            .id = core.ids.BondId.fromIndex(@intCast(bond_count)),
            .input_index = @intCast(bond_count),
            .start = core.ids.AtomId.fromIndex(start.index()),
            .end = core.ids.AtomId.fromIndex(end.index()),
            .input_order = .single,
            .effective_order = .single,
        };
        bond_count += 1;
    }
    if (bond_count != unique_edge_count) return error.InvalidMapping;

    var meta_graph = try topology.Graph.init(allocator, meta_atoms, meta_bonds);
    defer meta_graph.deinit();
    if (meta_graph.component_count != 1) return null;
    var meta_rings = try topology.RingMembership.init(allocator, meta_graph, meta_bonds);
    defer meta_rings.deinit();
    var fragmentation = try layout.Fragmentation.init(allocator, meta_atoms, meta_bonds, meta_graph, meta_rings);
    defer fragmentation.deinit();
    try layout.initializeCoordinates(allocator, meta_atoms, meta_bonds, meta_graph, meta_rings, fragmentation);
    try orientComponents(allocator, meta_atoms, meta_bonds, meta_graph, meta_rings, null);

    const coordinates = allocator.alloc(Vec2, meta_atoms.len) catch return error.OutOfMemory;
    for (meta_atoms, coordinates) |atom, *coordinate| coordinate.* = atom.coordinates;
    return coordinates;
}

fn rotateGeneralComponent(
    allocator: std.mem.Allocator,
    atoms: []model.Atom,
    graph: topology.Graph,
    rings: topology.RingMembership,
    relations: []const ProximityRelation,
    templates: []const Vec2,
    component: core.ids.MoleculeId,
) core.errors.Error!bool {
    const members = graph.componentMembers(component);
    if (members.len < 2) return true;
    const neighbors = allocator.alloc(core.ids.MoleculeId, graph.component_count) catch return error.OutOfMemory;
    defer allocator.free(neighbors);
    var neighbor_count: usize = 0;
    for (relations) |relation| {
        const local = localEndpoint(graph, relation, component) orelse continue;
        const other_atom = if (local == relation.start) relation.end else relation.start;
        const other = graph.component(other_atom) orelse return error.InvalidMapping;
        if (other == component or std.mem.indexOfScalar(core.ids.MoleculeId, neighbors[0..neighbor_count], other) != null) continue;
        neighbors[neighbor_count] = other;
        neighbor_count += 1;
    }
    if (neighbor_count == 0) return false;
    if (neighbor_count == 1) {
        const site = try proximitySiteToward(atoms, graph, rings, relations, component, neighbors[0]);
        if (site.count == 0) return false;
        const direction = geometry.subtract(templates[component.index()], templates[neighbors[0].index()]);
        const angle = geometry.signedAngle(geometry.scale(site.addition, -1), .{}, direction) * (-std.math.pi / 180.0);
        const sine = std.math.sin(angle);
        const cosine = std.math.cos(angle);
        for (members) |atom| {
            const relative = geometry.subtract(atoms[atom.index()].coordinates, site.center);
            atoms[atom.index()].coordinates = geometry.add(site.center, geometry.rotate(relative, sine, cosine));
        }
        return true;
    }

    const template_vectors = allocator.alloc(Vec2, neighbor_count) catch return error.OutOfMemory;
    defer allocator.free(template_vectors);
    const addition_vectors = allocator.alloc(Vec2, neighbor_count) catch return error.OutOfMemory;
    defer allocator.free(addition_vectors);
    for (neighbors[0..neighbor_count], template_vectors, addition_vectors) |neighbor, *template, *addition| {
        template.* = geometry.subtract(templates[neighbor.index()], templates[component.index()]);
        const site = try proximitySiteToward(atoms, graph, rings, relations, component, neighbor);
        if (site.count == 0) return false;
        addition.* = site.addition;
    }
    const matrix = try layout.templates.alignmentMatrix(template_vectors, addition_vectors);
    var center: Vec2 = .{};
    for (members) |atom| center = geometry.add(center, atoms[atom.index()].coordinates);
    center = geometry.scale(center, 1 / @as(f32, @floatFromInt(members.len)));
    for (members) |atom| {
        const relative = geometry.subtract(atoms[atom.index()].coordinates, center);
        atoms[atom.index()].coordinates = geometry.add(center, .{
            .x = relative.x * matrix[0] + relative.y * matrix[1],
            .y = relative.x * matrix[2] + relative.y * matrix[3],
        });
    }
    return true;
}

const ProximitySite = struct {
    center: Vec2,
    addition: Vec2,
    count: u32,
};

fn proximitySiteToward(
    atoms: []const model.Atom,
    graph: topology.Graph,
    rings: topology.RingMembership,
    relations: []const ProximityRelation,
    component: core.ids.MoleculeId,
    neighbor: core.ids.MoleculeId,
) core.errors.Error!ProximitySite {
    var result = ProximitySite{ .center = .{}, .addition = .{}, .count = 0 };
    for (relations) |relation| {
        const local = localEndpoint(graph, relation, component) orelse continue;
        const other = if (local == relation.start) relation.end else relation.start;
        if (graph.component(other) != neighbor) continue;
        result.center = geometry.add(result.center, atoms[local.index()].coordinates);
        result.addition = geometry.add(result.addition, singleAdditionVector(atoms, graph, rings, local));
        result.count += 1;
    }
    if (result.count != 0) result.center = geometry.scale(result.center, 1 / @as(f32, @floatFromInt(result.count)));
    result.addition = geometry.normalize(result.addition);
    return result;
}

fn proximitySite(
    atoms: []const model.Atom,
    graph: topology.Graph,
    rings: topology.RingMembership,
    relations: []const ProximityRelation,
    component: core.ids.MoleculeId,
) core.errors.Error!ProximitySite {
    var result = ProximitySite{ .center = .{}, .addition = .{}, .count = 0 };
    for (relations) |relation| {
        const local = localEndpoint(graph, relation, component) orelse continue;
        const other = if (local == relation.start) relation.end else relation.start;
        const other_component = graph.component(other) orelse return error.InvalidMapping;
        if (other_component == component) continue;
        result.center = geometry.add(result.center, atoms[local.index()].coordinates);
        result.addition = geometry.add(result.addition, singleAdditionVector(atoms, graph, rings, local));
        result.count += 1;
    }
    if (result.count != 0) result.center = geometry.scale(result.center, 1 / @as(f32, @floatFromInt(result.count)));
    result.addition = geometry.normalize(result.addition);
    return result;
}

fn replaceTerminalSingletons(
    atoms: []model.Atom,
    graph: topology.Graph,
    rings: topology.RingMembership,
    relations: []const ProximityRelation,
    counter: u32,
) core.errors.Error!void {
    for (0..graph.component_count) |raw_component| {
        const component = core.ids.MoleculeId.fromIndex(@intCast(raw_component));
        const members = graph.componentMembers(component);
        if (members.len != 1 or try componentNeighborCount(graph, relations, component) != 1) continue;
        var coordinates: Vec2 = .{};
        var count: u32 = 0;
        var fallback: ?Vec2 = null;
        for (relations) |relation| {
            const local = localEndpoint(graph, relation, component) orelse continue;
            const partner = if (local == relation.start) relation.end else relation.start;
            if (graph.component(partner) == component) continue;
            if (fallback == null) fallback = atoms[partner.index()].coordinates;
            var addition = singleAdditionVector(atoms, graph, rings, partner);
            if (geometry.length(addition) < geometry.epsilon) continue;
            addition = geometry.scale(geometry.normalize(addition), core.math.bond_length * @as(f32, @floatFromInt(counter)));
            coordinates = geometry.add(coordinates, geometry.add(atoms[partner.index()].coordinates, addition));
            count += 1;
        }
        atoms[members[0].index()].coordinates = if (count != 0)
            geometry.scale(coordinates, 1 / @as(f32, @floatFromInt(count)))
        else
            fallback orelse atoms[members[0].index()].coordinates;
    }
}

fn componentNeighborCount(
    graph: topology.Graph,
    relations: []const ProximityRelation,
    component: core.ids.MoleculeId,
) core.errors.Error!usize {
    var count: usize = 0;
    for (0..graph.component_count) |raw_other| {
        const other = core.ids.MoleculeId.fromIndex(@intCast(raw_other));
        if (other == component) continue;
        for (relations) |relation| {
            const start = graph.component(relation.start) orelse return error.InvalidMapping;
            const end = graph.component(relation.end) orelse return error.InvalidMapping;
            if ((start == component and end == other) or (end == component and start == other)) {
                count += 1;
                break;
            }
        }
    }
    return count;
}

fn componentsClash(atoms: []const model.Atom, graph: topology.Graph, distance: f32) bool {
    const squared_distance = distance * distance;
    for (0..graph.component_count) |raw_first| {
        const first = graph.componentMembers(core.ids.MoleculeId.fromIndex(@intCast(raw_first)));
        for (raw_first + 1..graph.component_count) |raw_second| {
            const second = graph.componentMembers(core.ids.MoleculeId.fromIndex(@intCast(raw_second)));
            for (first) |first_atom| {
                for (second) |second_atom| {
                    if (geometry.squaredLength(geometry.subtract(
                        atoms[first_atom.index()].coordinates,
                        atoms[second_atom.index()].coordinates,
                    )) < squared_distance) return true;
                }
            }
        }
    }
    return false;
}

fn singleAdditionVector(
    atoms: []const model.Atom,
    graph: topology.Graph,
    rings: topology.RingMembership,
    atom: core.ids.AtomId,
) Vec2 {
    var output: Vec2 = .{};
    var total_weight: f32 = 0;
    for (graph.neighbors(atom)) |neighbor| {
        const weight: f32 = if (shareRing(rings, atom, neighbor)) 4 else 1;
        output = geometry.add(output, geometry.scale(
            geometry.subtract(atoms[neighbor.index()].coordinates, atoms[atom.index()].coordinates),
            weight,
        ));
        total_weight += weight;
    }
    if (total_weight != 0) output = geometry.scale(output, 1 / total_weight);
    return geometry.scale(output, -1);
}

fn shareRing(rings: topology.RingMembership, first: core.ids.AtomId, second: core.ids.AtomId) bool {
    for (rings.atomRings(first)) |first_ring| {
        for (rings.atomRings(second)) |second_ring| if (first_ring == second_ring) return true;
    }
    return false;
}

/// Place the largest component at its existing position, then place neutral
/// and charged components in upstream molecule order on the first
/// non-clashing point of the pinned ten-level, one-bond-length grid.
/// Proximity-related placement is owned by a later controller slice.
pub fn arrangeComponents(
    allocator: std.mem.Allocator,
    atoms: []model.Atom,
    graph: topology.Graph,
) core.errors.Error!void {
    return arrangeComponentsExcluding(allocator, atoms, graph, &.{});
}

/// Arrange ordinary components while leaving residue-only components for the
/// dedicated crown/protein placement phase. A mixed excluded/non-excluded
/// component violates the residue representative ownership contract.
pub fn arrangeComponentsExcluding(
    allocator: std.mem.Allocator,
    atoms: []model.Atom,
    graph: topology.Graph,
    excluded_atoms: []const bool,
) core.errors.Error!void {
    if (excluded_atoms.len != 0 and excluded_atoms.len != atoms.len) return error.InvalidMapping;
    var excluded_components: []bool = &.{};
    defer if (excluded_components.len != 0) allocator.free(excluded_components);
    var active_count: usize = graph.component_count;
    if (excluded_atoms.len != 0) {
        excluded_components = allocator.alloc(bool, graph.component_count) catch return error.OutOfMemory;
        active_count = 0;
        for (excluded_components, 0..) |*excluded, raw_component| {
            const members = graph.componentMembers(core.ids.MoleculeId.fromIndex(@intCast(raw_component)));
            var excluded_count: usize = 0;
            for (members) |atom| excluded_count += @intFromBool(excluded_atoms[atom.index()]);
            if (excluded_count != 0 and excluded_count != members.len) return error.InvalidMapping;
            excluded.* = excluded_count == members.len;
            active_count += @intFromBool(!excluded.*);
        }
    }
    if (active_count < 2) return;
    const placed = allocator.alloc(bool, graph.component_count) catch return error.OutOfMemory;
    defer allocator.free(placed);
    @memset(placed, false);
    const charged_atom_used = allocator.alloc(bool, atoms.len) catch return error.OutOfMemory;
    defer allocator.free(charged_atom_used);
    @memset(charged_atom_used, false);

    var central_index: u32 = if (excluded_components.len == 0)
        0
    else
        @intCast(std.mem.indexOfScalar(bool, excluded_components, false).?);
    var central_size = graph.componentMembers(core.ids.MoleculeId.fromIndex(central_index)).len;
    for (0..graph.component_count) |raw_index| {
        if (excluded_components.len != 0 and excluded_components[raw_index]) continue;
        const index: u32 = @intCast(raw_index);
        const size = graph.componentMembers(core.ids.MoleculeId.fromIndex(index)).len;
        if (size > central_size) {
            central_index = index;
            central_size = size;
        }
    }
    placed[central_index] = true;
    const center = componentCenter(atoms, graph.componentMembers(core.ids.MoleculeId.fromIndex(central_index)));

    for (0..graph.component_count) |raw_index| {
        const index: u32 = @intCast(raw_index);
        if ((excluded_components.len != 0 and excluded_components[index]) or placed[index]) continue;
        const component = core.ids.MoleculeId.fromIndex(index);
        const members = graph.componentMembers(component);
        var charge: i64 = 0;
        for (members) |atom| charge += atoms[atom.index()].formal_charge;
        if (charge != 0) continue;

        const component_bounds = componentBounds(atoms, members);
        const half_width = (component_bounds.max.x - component_bounds.min.x) * 0.5;
        const half_height = (component_bounds.max.y - component_bounds.min.y) * 0.5;
        const old_center = Vec2{
            .x = (component_bounds.max.x + component_bounds.min.x) * 0.5,
            .y = (component_bounds.max.y + component_bounds.min.y) * 0.5,
        };
        const target = findGridPoint(atoms, graph, placed, center, half_width, half_height, core.math.bond_length * 1.8);
        const translation = Vec2{ .x = target.x - old_center.x, .y = target.y - old_center.y };
        for (members) |atom| {
            atoms[atom.index()].coordinates.x += translation.x;
            atoms[atom.index()].coordinates.y += translation.y;
        }
        placed[index] = true;
    }

    for (0..graph.component_count) |raw_index| {
        const index: u32 = @intCast(raw_index);
        if ((excluded_components.len != 0 and excluded_components[index]) or placed[index]) continue;
        const component = core.ids.MoleculeId.fromIndex(index);
        const members = graph.componentMembers(component);
        var charge: i64 = 0;
        for (members) |atom| charge += atoms[atom.index()].formal_charge;
        if (charge == 0) continue;

        var search_center = center;
        outer: for (0..graph.component_count) |placed_index| {
            if (!placed[placed_index]) continue;
            for (graph.componentMembers(core.ids.MoleculeId.fromIndex(@intCast(placed_index)))) |atom| {
                const atom_index = atom.index();
                const atom_charge = atoms[atom_index].formal_charge;
                if (atom_charge == 0 or charged_atom_used[atom_index]) continue;
                if (@as(i64, atom_charge) * charge < 0) {
                    charged_atom_used[atom_index] = true;
                    search_center = atoms[atom_index].coordinates;
                    break :outer;
                }
            }
        }

        const component_bounds = componentBounds(atoms, members);
        const half_width = (component_bounds.max.x - component_bounds.min.x) * 0.5;
        const half_height = (component_bounds.max.y - component_bounds.min.y) * 0.5;
        const old_center = Vec2{
            .x = (component_bounds.max.x + component_bounds.min.x) * 0.5,
            .y = (component_bounds.max.y + component_bounds.min.y) * 0.5,
        };
        const target = findGridPoint(atoms, graph, placed, search_center, half_width, half_height, core.math.bond_length * 0.8);
        const translation = Vec2{ .x = target.x - old_center.x, .y = target.y - old_center.y };
        for (members) |atom| {
            atoms[atom.index()].coordinates.x += translation.x;
            atoms[atom.index()].coordinates.y += translation.y;
        }
        placed[index] = true;
    }
}

fn findGridPoint(
    atoms: []const model.Atom,
    graph: topology.Graph,
    placed: []const bool,
    center: Vec2,
    half_width: f32,
    half_height: f32,
    clearance: f32,
) Vec2 {
    const step = core.math.bond_length;
    for (0..10) |level| {
        const distance = @as(f32, @floatFromInt(level + 1)) * step;
        const right = Vec2{ .x = center.x + distance, .y = center.y };
        const left = Vec2{ .x = center.x - distance, .y = center.y };
        const bottom = Vec2{ .x = center.x, .y = center.y - distance };
        const top = Vec2{ .x = center.x, .y = center.y + distance };
        const cardinal = [_]Vec2{ center, right, left, bottom, top };
        for (cardinal) |candidate| if (pointIsClear(atoms, graph, placed, candidate, half_width, half_height, clearance)) return candidate;
        for (0..level) |offset_index| {
            const offset = @as(f32, @floatFromInt(offset_index + 1)) * step;
            const edge = [_]Vec2{
                .{ .x = right.x, .y = right.y + offset },
                .{ .x = right.x, .y = right.y - offset },
                .{ .x = left.x, .y = left.y + offset },
                .{ .x = left.x, .y = left.y - offset },
                .{ .x = bottom.x + offset, .y = bottom.y },
                .{ .x = bottom.x - offset, .y = bottom.y },
                .{ .x = top.x + offset, .y = top.y },
                .{ .x = top.x - offset, .y = top.y },
            };
            for (edge) |candidate| if (pointIsClear(atoms, graph, placed, candidate, half_width, half_height, clearance)) return candidate;
        }
        const corners = [_]Vec2{
            .{ .x = center.x + distance, .y = center.y + distance },
            .{ .x = center.x + distance, .y = center.y - distance },
            .{ .x = center.x - distance, .y = center.y + distance },
            .{ .x = center.x - distance, .y = center.y - distance },
        };
        for (corners) |candidate| if (pointIsClear(atoms, graph, placed, candidate, half_width, half_height, clearance)) return candidate;
    }
    return center;
}

fn pointIsClear(
    atoms: []const model.Atom,
    graph: topology.Graph,
    placed: []const bool,
    candidate: Vec2,
    half_width: f32,
    half_height: f32,
    clearance: f32,
) bool {
    for (0..graph.component_count) |raw_index| {
        if (!placed[raw_index]) continue;
        for (graph.componentMembers(core.ids.MoleculeId.fromIndex(@intCast(raw_index)))) |atom| {
            const coordinate = atoms[atom.index()].coordinates;
            if (coordinate.x < candidate.x + clearance + half_width and
                coordinate.x > candidate.x - clearance - half_width and
                coordinate.y < candidate.y + clearance + half_height and
                coordinate.y > candidate.y - clearance - half_height)
            {
                return false;
            }
        }
    }
    return true;
}

fn componentBounds(atoms: []const model.Atom, members: []const core.ids.AtomId) Bounds {
    var result = Bounds{ .min = atoms[members[0].index()].coordinates, .max = atoms[members[0].index()].coordinates };
    for (members) |atom| {
        const coordinate = atoms[atom.index()].coordinates;
        if (coordinate.x < result.min.x) result.min.x = coordinate.x;
        if (coordinate.y < result.min.y) result.min.y = coordinate.y;
        if (coordinate.x > result.max.x) result.max.x = coordinate.x;
        if (coordinate.y > result.max.y) result.max.y = coordinate.y;
    }
    return result;
}

fn componentCenter(atoms: []const model.Atom, members: []const core.ids.AtomId) Vec2 {
    var result: Vec2 = .{};
    for (members) |atom| {
        result.x += atoms[atom.index()].coordinates.x;
        result.y += atoms[atom.index()].coordinates.y;
    }
    result.x /= @floatFromInt(members.len);
    result.y /= @floatFromInt(members.len);
    return result;
}

fn testAtom(index: u32, x: f32) model.Atom {
    return .{
        .id = core.ids.AtomId.fromIndex(index),
        .input_index = index,
        .atomic_number = .carbon,
        .coordinates = .{ .x = x },
    };
}

test "neutral components use largest-first center and pinned grid order" {
    var atoms = [_]model.Atom{ testAtom(0, 0), testAtom(1, 50), testAtom(2, 0) };
    const bonds = [_]model.Bond{.{
        .id = core.ids.BondId.fromIndex(0),
        .input_index = 0,
        .start = core.ids.AtomId.fromIndex(0),
        .end = core.ids.AtomId.fromIndex(1),
        .input_order = .single,
        .effective_order = .single,
    }};
    var graph = try topology.Graph.init(std.testing.allocator, &atoms, &bonds);
    defer graph.deinit();
    try arrangeComponents(std.testing.allocator, &atoms, graph);
    try std.testing.expectEqual(Vec2{ .x = 0, .y = 0 }, atoms[0].coordinates);
    try std.testing.expectEqual(Vec2{ .x = 50, .y = 0 }, atoms[1].coordinates);
    // The first level clashes throughout. At the second level +x and -x
    // still clash; -y is the first clear point in pinned candidate order.
    try std.testing.expectEqual(Vec2{ .x = 25, .y = -100 }, atoms[2].coordinates);
}

test "counterion uses first unused opposite charged atom and shorter clearance" {
    var atoms = [_]model.Atom{ testAtom(0, 0), testAtom(1, 50), testAtom(2, 0) };
    atoms[0].formal_charge = 1;
    atoms[2].formal_charge = -1;
    const bonds = [_]model.Bond{.{
        .id = core.ids.BondId.fromIndex(0),
        .input_index = 0,
        .start = core.ids.AtomId.fromIndex(0),
        .end = core.ids.AtomId.fromIndex(1),
        .input_order = .single,
        .effective_order = .single,
    }};
    var graph = try topology.Graph.init(std.testing.allocator, &atoms, &bonds);
    defer graph.deinit();
    try arrangeComponents(std.testing.allocator, &atoms, graph);
    try std.testing.expectEqual(Vec2{ .x = -50, .y = 0 }, atoms[2].coordinates);
}

test "proximity center uses relation count then size with first-wins ties" {
    var atoms = [_]model.Atom{
        testAtom(0, 0), testAtom(1, 50), testAtom(2, 0), testAtom(3, 50), testAtom(4, 0),
    };
    const bonds = [_]model.Bond{
        .{
            .id = core.ids.BondId.fromIndex(0),
            .input_index = 0,
            .start = core.ids.AtomId.fromIndex(0),
            .end = core.ids.AtomId.fromIndex(1),
            .input_order = .single,
            .effective_order = .single,
        },
        .{
            .id = core.ids.BondId.fromIndex(1),
            .input_index = 1,
            .start = core.ids.AtomId.fromIndex(2),
            .end = core.ids.AtomId.fromIndex(3),
            .input_order = .single,
            .effective_order = .single,
        },
    };
    var graph = try topology.Graph.init(std.testing.allocator, &atoms, &bonds);
    defer graph.deinit();
    const first_relation = [_]ProximityRelation{.{
        .start = core.ids.AtomId.fromIndex(0),
        .end = core.ids.AtomId.fromIndex(4),
    }};
    try std.testing.expectEqual(core.ids.MoleculeId.fromIndex(0), try selectProximityCenter(graph, &first_relation));
    const more_relations = [_]ProximityRelation{
        .{ .start = core.ids.AtomId.fromIndex(2), .end = core.ids.AtomId.fromIndex(4) },
        .{ .start = core.ids.AtomId.fromIndex(3), .end = core.ids.AtomId.fromIndex(4) },
        .{ .start = core.ids.AtomId.fromIndex(2), .end = core.ids.AtomId.fromIndex(0) },
    };
    try std.testing.expectEqual(core.ids.MoleculeId.fromIndex(1), try selectProximityCenter(graph, &more_relations));
}

test "residue interactions are reserved for residue placement" {
    const residues = [_]model.Residue{.{
        .id = core.ids.ResidueId.fromIndex(0),
        .atom = core.ids.AtomId.fromIndex(2),
        .chain_start = 0,
        .chain_len = 1,
        .residue_number = 7,
    }};
    const endpoints = ProximityRelation{
        .start = core.ids.AtomId.fromIndex(1),
        .end = core.ids.AtomId.fromIndex(2),
    };
    try std.testing.expect(isResidueRelation(&residues, .{ .residue_interaction = 0 }, endpoints));
    try std.testing.expect(!isResidueRelation(&residues, .{ .bond = core.ids.BondId.fromIndex(0) }, endpoints));
    try std.testing.expect(!isResidueRelation(&.{}, .{ .residue_interaction = 0 }, endpoints));
}

test "isolated acyclic bond rotates to the upstream horizontal preference" {
    var atoms = [_]model.Atom{ testAtom(0, 0), testAtom(1, 0) };
    atoms[1].coordinates = .{ .y = 50 };
    const bonds = [_]model.Bond{.{
        .id = core.ids.BondId.fromIndex(0),
        .input_index = 0,
        .start = core.ids.AtomId.fromIndex(0),
        .end = core.ids.AtomId.fromIndex(1),
        .input_order = .single,
        .effective_order = .single,
    }};
    var graph = try topology.Graph.init(std.testing.allocator, &atoms, &bonds);
    defer graph.deinit();
    var rings = try topology.RingMembership.init(std.testing.allocator, graph, &bonds);
    defer rings.deinit();
    try orientComponents(std.testing.allocator, &atoms, &bonds, graph, rings, null);
    try std.testing.expectApproxEqAbs(atoms[0].coordinates.y, atoms[1].coordinates.y, 0.05);
    try std.testing.expectApproxEqAbs(core.math.bond_length, @abs(atoms[1].coordinates.x - atoms[0].coordinates.x), 0.001);
}

test "first crossing proximity pair mirrors only its local component" {
    var atoms = [_]model.Atom{
        testAtom(0, -1), testAtom(1, 1), testAtom(2, 1), testAtom(3, -1),
    };
    atoms[2].coordinates.y = 2;
    atoms[3].coordinates.y = 2;
    const bonds = [_]model.Bond{.{
        .id = core.ids.BondId.fromIndex(0),
        .input_index = 0,
        .start = core.ids.AtomId.fromIndex(0),
        .end = core.ids.AtomId.fromIndex(1),
        .input_order = .single,
        .effective_order = .single,
    }};
    var graph = try topology.Graph.init(std.testing.allocator, &atoms, &bonds);
    defer graph.deinit();
    const relations = [_]ProximityRelation{
        .{ .start = core.ids.AtomId.fromIndex(0), .end = core.ids.AtomId.fromIndex(2) },
        .{ .start = core.ids.AtomId.fromIndex(1), .end = core.ids.AtomId.fromIndex(3) },
    };
    const placed = [_]bool{ false, true, true };
    try std.testing.expect(try flipFirstCrossingRelations(&atoms, graph, &relations, core.ids.MoleculeId.fromIndex(0), &placed));
    try std.testing.expectEqual(Vec2{ .x = 1 }, atoms[0].coordinates);
    try std.testing.expectEqual(Vec2{ .x = -1 }, atoms[1].coordinates);
    try std.testing.expectEqual(Vec2{ .x = 1, .y = 2 }, atoms[2].coordinates);
    try std.testing.expectEqual(Vec2{ .x = -1, .y = 2 }, atoms[3].coordinates);
}

test "proximity child aligns interaction site with parent free valence" {
    var atoms = [_]model.Atom{ testAtom(0, 0), testAtom(1, 50), testAtom(2, 0), testAtom(3, 50) };
    const bonds = [_]model.Bond{
        .{
            .id = core.ids.BondId.fromIndex(0),
            .input_index = 0,
            .start = core.ids.AtomId.fromIndex(0),
            .end = core.ids.AtomId.fromIndex(1),
            .input_order = .single,
            .effective_order = .single,
        },
        .{
            .id = core.ids.BondId.fromIndex(1),
            .input_index = 1,
            .start = core.ids.AtomId.fromIndex(2),
            .end = core.ids.AtomId.fromIndex(3),
            .input_order = .single,
            .effective_order = .single,
        },
    };
    var graph = try topology.Graph.init(std.testing.allocator, &atoms, &bonds);
    defer graph.deinit();
    var rings = try topology.RingMembership.init(std.testing.allocator, graph, &bonds);
    defer rings.deinit();
    const relations = [_]ProximityRelation{.{
        .start = core.ids.AtomId.fromIndex(1),
        .end = core.ids.AtomId.fromIndex(2),
    }};
    var placed = [_]bool{ true, false };
    try std.testing.expect(try placeProximityChild(
        &atoms,
        graph,
        rings,
        &relations,
        core.ids.MoleculeId.fromIndex(1),
        core.ids.MoleculeId.fromIndex(0),
        &placed,
    ));
    try std.testing.expectEqual(Vec2{ .x = 200 }, atoms[2].coordinates);
    try std.testing.expectEqual(Vec2{ .x = 250 }, atoms[3].coordinates);
    try std.testing.expect(placed[1]);
}

fn arrangeAndDiscard(allocator: std.mem.Allocator, atoms: []model.Atom, graph: topology.Graph) !void {
    try arrangeComponents(allocator, atoms, graph);
}

fn arrangeProximityExcludingAndDiscard(
    allocator: std.mem.Allocator,
    source_atoms: []const model.Atom,
    bonds: []const model.Bond,
    relations: []const ProximityRelation,
    excluded: []const bool,
) !void {
    const atoms = try allocator.dupe(model.Atom, source_atoms);
    defer allocator.free(atoms);
    if (!try arrangeProximityComponentsExcluding(allocator, atoms, bonds, relations, excluded)) return error.InvalidMapping;
}

test "neutral component placement cleans every allocation failure" {
    var atoms = [_]model.Atom{ testAtom(0, 0), testAtom(1, 50), testAtom(2, 0) };
    const bonds = [_]model.Bond{.{
        .id = core.ids.BondId.fromIndex(0),
        .input_index = 0,
        .start = core.ids.AtomId.fromIndex(0),
        .end = core.ids.AtomId.fromIndex(1),
        .input_order = .single,
        .effective_order = .single,
    }};
    var graph = try topology.Graph.init(std.testing.allocator, &atoms, &bonds);
    defer graph.deinit();
    try core.oom.checkAllocationFailures(std.testing.allocator, arrangeAndDiscard, .{ &atoms, graph });
}

test "component arrangement ignores residue-only components" {
    var atoms = [_]model.Atom{ testAtom(0, 0), testAtom(1, 50), testAtom(2, 0), testAtom(3, 17) };
    atoms[3].coordinates.y = 23;
    const bonds = [_]model.Bond{.{
        .id = core.ids.BondId.fromIndex(0),
        .input_index = 0,
        .start = core.ids.AtomId.fromIndex(0),
        .end = core.ids.AtomId.fromIndex(1),
        .input_order = .single,
        .effective_order = .single,
    }};
    var graph = try topology.Graph.init(std.testing.allocator, &atoms, &bonds);
    defer graph.deinit();
    const excluded = [_]bool{ false, false, false, true };
    try arrangeComponentsExcluding(std.testing.allocator, &atoms, graph, &excluded);
    try std.testing.expect(atoms[2].coordinates.x != 0 or atoms[2].coordinates.y != 0);
    try std.testing.expectEqual(Vec2{ .x = 17, .y = 23 }, atoms[3].coordinates);

    const mixed = [_]bool{ true, false, false, false };
    try std.testing.expectError(error.InvalidMapping, arrangeComponentsExcluding(std.testing.allocator, &atoms, graph, &mixed));
}

test "proximity arrangement compacts excluded residue components and cleans allocation failures" {
    const source_atoms = [_]model.Atom{ testAtom(0, 0), testAtom(1, 50), testAtom(2, 0), testAtom(3, 19) };
    const bonds = [_]model.Bond{
        .{
            .id = core.ids.BondId.fromIndex(0),
            .input_index = 0,
            .start = core.ids.AtomId.fromIndex(0),
            .end = core.ids.AtomId.fromIndex(1),
            .input_order = .single,
            .effective_order = .single,
        },
        .{
            .id = core.ids.BondId.fromIndex(1),
            .input_index = 1,
            .start = core.ids.AtomId.fromIndex(1),
            .end = core.ids.AtomId.fromIndex(2),
            .input_order = .zero,
            .effective_order = .zero,
        },
    };
    const relations = [_]ProximityRelation{.{
        .start = core.ids.AtomId.fromIndex(1),
        .end = core.ids.AtomId.fromIndex(2),
    }};
    const excluded = [_]bool{ false, false, false, true };
    try core.oom.checkAllocationFailures(std.testing.allocator, arrangeProximityExcludingAndDiscard, .{
        &source_atoms,
        &bonds,
        &relations,
        &excluded,
    });
}
