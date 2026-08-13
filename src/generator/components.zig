const std = @import("std");
const core = @import("core");
const model = @import("model");
const geometry = @import("geometry");
const topology = @import("topology");

const Vec2 = core.math.Vec2;
const Bounds = struct { min: Vec2, max: Vec2 };

pub const ProximityRelation = struct {
    start: core.ids.AtomId,
    end: core.ids.AtomId,
};

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

/// Apply upstream's global bestRotation/maybeFlip rules to acyclic components.
/// Ring-fusion and peptide candidate bonuses are deferred until those
/// controller records are available; affected components are left untouched.
pub fn orientAcyclicComponents(
    allocator: std.mem.Allocator,
    atoms: []model.Atom,
    bonds: []const model.Bond,
    graph: topology.Graph,
    rings: topology.RingMembership,
) core.errors.Error!void {
    for (0..graph.component_count) |raw_component| {
        const component = core.ids.MoleculeId.fromIndex(@intCast(raw_component));
        const members = graph.componentMembers(component);
        var has_ring = false;
        var constrained = false;
        for (members) |atom| {
            has_ring = has_ring or rings.atomRings(atom).len != 0;
            constrained = constrained or atoms[atom.index()].fixed or atoms[atom.index()].constrained;
        }
        if (has_ring or constrained) continue;

        var angles: std.ArrayList(WeightedAngle) = .empty;
        defer angles.deinit(allocator);
        for (members) |atom| {
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
        maybeFlipComponent(atoms, bonds, graph, component, members);
    }
}

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

fn maybeFlipComponent(atoms: []model.Atom, bonds: []const model.Bond, graph: topology.Graph, component: core.ids.MoleculeId, members: []const core.ids.AtomId) void {
    if (members.len < 2) return;
    const center = componentCenter(atoms, members);
    const bounds = componentBounds(atoms, members);
    const midpoint = Vec2{ .x = (bounds.max.x + bounds.min.x) * 0.5, .y = (bounds.max.y + bounds.min.y) * 0.5 };
    var score_x: f32 = 0;
    var score_y: f32 = 0;
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

/// Place the largest component at its existing position, then place neutral
/// and charged components in upstream molecule order on the first
/// non-clashing point of the pinned ten-level, one-bond-length grid.
/// Proximity-related placement is owned by a later controller slice.
pub fn arrangeComponents(
    allocator: std.mem.Allocator,
    atoms: []model.Atom,
    graph: topology.Graph,
) core.errors.Error!void {
    if (graph.component_count < 2) return;
    const placed = allocator.alloc(bool, graph.component_count) catch return error.OutOfMemory;
    defer allocator.free(placed);
    @memset(placed, false);
    const charged_atom_used = allocator.alloc(bool, atoms.len) catch return error.OutOfMemory;
    defer allocator.free(charged_atom_used);
    @memset(charged_atom_used, false);

    var central_index: u32 = 0;
    var central_size = graph.componentMembers(core.ids.MoleculeId.fromIndex(0)).len;
    for (1..graph.component_count) |raw_index| {
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
        if (placed[index]) continue;
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
        if (placed[index]) continue;
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
    try orientAcyclicComponents(std.testing.allocator, &atoms, &bonds, graph, rings);
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

fn arrangeAndDiscard(allocator: std.mem.Allocator, atoms: []model.Atom, graph: topology.Graph) !void {
    try arrangeComponents(allocator, atoms, graph);
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
    try std.testing.checkAllAllocationFailures(std.testing.allocator, arrangeAndDiscard, .{ &atoms, graph });
}
