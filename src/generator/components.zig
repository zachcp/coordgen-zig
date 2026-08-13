const std = @import("std");
const core = @import("core");
const model = @import("model");
const topology = @import("topology");

const Vec2 = core.math.Vec2;
const Bounds = struct { min: Vec2, max: Vec2 };

/// Place the largest component at its existing position, then place neutral
/// components in upstream molecule order on the first non-clashing point of
/// the pinned ten-level, one-bond-length grid. Charged and proximity-related
/// placement is owned by later controller slices.
pub fn arrangeNeutralComponents(
    allocator: std.mem.Allocator,
    atoms: []model.Atom,
    graph: topology.Graph,
) core.errors.Error!void {
    if (graph.component_count < 2) return;
    const placed = allocator.alloc(bool, graph.component_count) catch return error.OutOfMemory;
    defer allocator.free(placed);
    @memset(placed, false);

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
        const target = findGridPoint(atoms, graph, placed, center, half_width, half_height);
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
) Vec2 {
    const step = core.math.bond_length;
    for (0..10) |level| {
        const distance = @as(f32, @floatFromInt(level + 1)) * step;
        const right = Vec2{ .x = center.x + distance, .y = center.y };
        const left = Vec2{ .x = center.x - distance, .y = center.y };
        const bottom = Vec2{ .x = center.x, .y = center.y - distance };
        const top = Vec2{ .x = center.x, .y = center.y + distance };
        const cardinal = [_]Vec2{ center, right, left, bottom, top };
        for (cardinal) |candidate| if (pointIsClear(atoms, graph, placed, candidate, half_width, half_height)) return candidate;
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
            for (edge) |candidate| if (pointIsClear(atoms, graph, placed, candidate, half_width, half_height)) return candidate;
        }
        const corners = [_]Vec2{
            .{ .x = center.x + distance, .y = center.y + distance },
            .{ .x = center.x + distance, .y = center.y - distance },
            .{ .x = center.x - distance, .y = center.y + distance },
            .{ .x = center.x - distance, .y = center.y - distance },
        };
        for (corners) |candidate| if (pointIsClear(atoms, graph, placed, candidate, half_width, half_height)) return candidate;
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
) bool {
    const clearance = core.math.bond_length * 1.8;
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
    try arrangeNeutralComponents(std.testing.allocator, &atoms, graph);
    try std.testing.expectEqual(Vec2{ .x = 0, .y = 0 }, atoms[0].coordinates);
    try std.testing.expectEqual(Vec2{ .x = 50, .y = 0 }, atoms[1].coordinates);
    // The first level clashes throughout. At the second level +x and -x
    // still clash; -y is the first clear point in pinned candidate order.
    try std.testing.expectEqual(Vec2{ .x = 25, .y = -100 }, atoms[2].coordinates);
}

fn arrangeAndDiscard(allocator: std.mem.Allocator, atoms: []model.Atom, graph: topology.Graph) !void {
    try arrangeNeutralComponents(allocator, atoms, graph);
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
