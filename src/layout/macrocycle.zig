const std = @import("std");
const core = @import("core");
const geometry = @import("geometry");

pub const max_macrocycles: usize = 40;
pub const path_failed: i32 = -1000;
pub const sqrt_three_halves: f32 = 0.8660254037844386;

pub const HexCoords = struct {
    x: i32,
    y: i32,

    pub fn z(self: HexCoords) i32 {
        return -self.x - self.y;
    }

    pub fn distanceFrom(self: HexCoords, origin: HexCoords) u32 {
        return @max(@abs(self.x - origin.x), @abs(self.y - origin.y), @abs(self.z() - origin.z()));
    }

    pub fn rotate(self: HexCoords) HexCoords {
        return .{ .x = self.x + self.y, .y = -self.x };
    }

    pub fn vertex(self: HexCoords) VertexCoords {
        return .{ .x = self.x, .y = self.y, .z = self.z() };
    }

    pub fn neighbors(self: HexCoords) [6]HexCoords {
        return .{
            .{ .x = self.x + 1, .y = self.y },
            .{ .x = self.x + 1, .y = self.y - 1 },
            .{ .x = self.x, .y = self.y - 1 },
            .{ .x = self.x - 1, .y = self.y },
            .{ .x = self.x - 1, .y = self.y + 1 },
            .{ .x = self.x, .y = self.y + 1 },
        };
    }
};

pub const VertexCoords = struct {
    x: i32,
    y: i32,
    z: i32,

    pub fn add(self: VertexCoords, other: VertexCoords) VertexCoords {
        return .{ .x = self.x + other.x, .y = self.y + other.y, .z = self.z + other.z };
    }

    pub fn subtract(self: VertexCoords, other: VertexCoords) VertexCoords {
        return .{ .x = self.x - other.x, .y = self.y - other.y, .z = self.z - other.z };
    }

    pub fn toCartesian(self: VertexCoords) core.math.Vec2 {
        return .{
            .x = core.math.bond_length * sqrt_three_halves * @as(f32, @floatFromInt(self.x - self.z)),
            .y = core.math.bond_length * @as(f32, @floatFromInt(-self.x + 2 * self.y - self.z)) * 0.5,
        };
    }
};

pub const Polyomino = struct {
    const IncidentHexes = struct { items: [3]usize, len: usize };

    allocator: std.mem.Allocator,
    hexes: std.ArrayList(HexCoords) = .empty,
    pentagon_vertices: std.ArrayList(VertexCoords) = .empty,

    pub fn init(allocator: std.mem.Allocator) Polyomino {
        return .{ .allocator = allocator };
    }

    pub fn deinit(self: *Polyomino) void {
        self.pentagon_vertices.deinit(self.allocator);
        self.hexes.deinit(self.allocator);
        self.* = undefined;
    }

    pub fn clone(self: Polyomino, allocator: std.mem.Allocator) core.errors.Error!Polyomino {
        var result = Polyomino.init(allocator);
        errdefer result.deinit();
        result.hexes.appendSlice(allocator, self.hexes.items) catch return error.OutOfMemory;
        result.pentagon_vertices.appendSlice(allocator, self.pentagon_vertices.items) catch return error.OutOfMemory;
        return result;
    }

    pub fn clear(self: *Polyomino) void {
        self.hexes.clearRetainingCapacity();
        self.pentagon_vertices.clearRetainingCapacity();
    }

    pub fn addHex(self: *Polyomino, coordinates: HexCoords) core.errors.Error!void {
        if (self.indexOf(coordinates) != null) return error.InvalidMapping;
        self.hexes.append(self.allocator, coordinates) catch return error.OutOfMemory;
    }

    pub fn removeHex(self: *Polyomino, coordinates: HexCoords) core.errors.Error!void {
        const index = self.indexOf(coordinates) orelse return error.InvalidMapping;
        _ = self.hexes.orderedRemove(index);
    }

    pub fn contains(self: Polyomino, coordinates: HexCoords) bool {
        return self.indexOf(coordinates) != null;
    }

    fn indexOf(self: Polyomino, coordinates: HexCoords) ?usize {
        for (self.hexes.items, 0..) |candidate, index| {
            if (std.meta.eql(candidate, coordinates)) return index;
        }
        return null;
    }

    pub fn countNeighbors(self: Polyomino, coordinates: HexCoords) usize {
        var count: usize = 0;
        for (coordinates.neighbors()) |neighbor| count += @intFromBool(self.contains(neighbor));
        return count;
    }

    fn incidentHexes(self: Polyomino, vertex: VertexCoords) core.errors.Error!IncidentHexes {
        const direction = vertex.x + vertex.y + vertex.z;
        if (direction != 1 and direction != -1) return error.InvalidCoordinate;
        const candidates = [_]HexCoords{
            .{ .x = vertex.x - direction, .y = vertex.y },
            .{ .x = vertex.x, .y = vertex.y - direction },
            .{ .x = vertex.x, .y = vertex.y },
        };
        var result = IncidentHexes{ .items = .{ 0, 0, 0 }, .len = 0 };
        for (candidates) |candidate| if (self.indexOf(candidate)) |index| {
            result.items[result.len] = index;
            result.len += 1;
        };
        return result;
    }

    pub fn hexagonsAtVertex(self: Polyomino, vertex: VertexCoords) core.errors.Error!usize {
        return (try self.incidentHexes(vertex)).len;
    }

    pub fn freeVertexNeighborPositions(
        self: Polyomino,
        vertex: VertexCoords,
        output: *[3]HexCoords,
    ) core.errors.Error!usize {
        const direction = vertex.x + vertex.y + vertex.z;
        if (direction != 1 and direction != -1) return error.InvalidCoordinate;
        const candidates = [_]HexCoords{
            .{ .x = vertex.x - direction, .y = vertex.y },
            .{ .x = vertex.x, .y = vertex.y - direction },
            .{ .x = vertex.x, .y = vertex.y },
        };
        var count: usize = 0;
        for (candidates) |candidate| {
            if (self.contains(candidate)) continue;
            output[count] = candidate;
            count += 1;
        }
        return count;
    }

    fn followingVertex(center: HexCoords, vertex: VertexCoords) core.errors.Error!VertexCoords {
        var dx = vertex.x - center.x;
        var dy = vertex.y - center.y;
        var dz = vertex.z - center.z();
        if (dx == 0 and dy == 0) {
            dx = -dz;
            dz = 0;
        } else if (dx == 0 and dz == 0) {
            dz = -dy;
            dy = 0;
        } else if (dy == 0 and dz == 0) {
            dy = -dx;
            dx = 0;
        } else return error.InvalidCoordinate;
        return .{ .x = center.x + dx, .y = center.y + dy, .z = center.z() + dz };
    }

    fn outerVertex(self: Polyomino) core.errors.Error!VertexCoords {
        for (self.hexes.items) |hex| {
            const vertex = VertexCoords{ .x = hex.x + 1, .y = hex.y, .z = hex.z() };
            if (try self.hexagonsAtVertex(vertex) == 1) return vertex;
        }
        return error.InvalidMapping;
    }

    pub fn path(self: Polyomino, allocator: std.mem.Allocator) core.errors.Error![]VertexCoords {
        if (self.hexes.items.len == 0) return allocator.alloc(VertexCoords, 0) catch return error.OutOfMemory;
        const first = try self.outerVertex();
        var current = first;
        const initial_neighbors = try self.incidentHexes(first);
        if (initial_neighbors.len != 1) return error.InvalidMapping;
        var current_hex = initial_neighbors.items[0];
        var next = try followingVertex(self.hexes.items[current_hex], current);
        var output: std.ArrayList(VertexCoords) = .empty;
        defer output.deinit(allocator);
        while (true) {
            if (!containsVertex(self.pentagon_vertices.items, current)) {
                output.append(allocator, current) catch return error.OutOfMemory;
            }
            current = next;
            const neighbors = try self.incidentHexes(current);
            if (neighbors.len == 0 or neighbors.len > 2) return error.InvalidMapping;
            if (neighbors.len == 2) current_hex = if (neighbors.items[0] != current_hex) neighbors.items[0] else neighbors.items[1];
            next = try followingVertex(self.hexes.items[current_hex], current);
            if (std.meta.eql(current, first)) break;
        }
        return output.toOwnedSlice(allocator) catch return error.OutOfMemory;
    }

    pub fn coordinatesOfSubstituent(self: Polyomino, position: VertexCoords) core.errors.Error!VertexCoords {
        const neighbors = try self.incidentHexes(position);
        if (neighbors.len == 0 or neighbors.len > 2) return error.InvalidMapping;
        if (neighbors.len == 1) {
            const parent = self.hexes.items[neighbors.items[0]].vertex();
            var difference = position.subtract(parent);
            const direction: i32 = if (difference.x + difference.y + difference.z > 0) 1 else -1;
            if (difference.x == 0) difference.x -= direction;
            if (difference.y == 0) difference.y -= direction;
            if (difference.z == 0) difference.z -= direction;
            return parent.add(difference);
        }
        const first = self.hexes.items[neighbors.items[0]].vertex();
        const second = self.hexes.items[neighbors.items[1]].vertex();
        return second.subtract(position.subtract(first));
    }

    pub fn isTheSameAs(self: Polyomino, other: Polyomino, allocator: std.mem.Allocator) core.errors.Error!bool {
        if (self.hexes.items.len != other.hexes.items.len) return false;
        if (self.hexes.items.len == 0) return true;
        const target = allocator.dupe(HexCoords, other.hexes.items) catch return error.OutOfMemory;
        defer allocator.free(target);
        var lowest_x = self.hexes.items[0].x;
        var lowest_y = self.hexes.items[0].y;
        for (self.hexes.items) |hex| {
            lowest_x = @min(lowest_x, hex.x);
            lowest_y = @min(lowest_y, hex.y);
        }
        for (0..6) |_| {
            var target_x = target[0].x;
            var target_y = target[0].y;
            for (target) |hex| {
                target_x = @min(target_x, hex.x);
                target_y = @min(target_y, hex.y);
            }
            var same = true;
            for (target) |*hex| {
                hex.x += lowest_x - target_x;
                hex.y += lowest_y - target_y;
                if (!self.contains(hex.*)) same = false;
            }
            if (same) return true;
            for (target) |*hex| hex.* = hex.rotate();
        }
        return false;
    }

    pub fn allFreeNeighbors(self: Polyomino, allocator: std.mem.Allocator) core.errors.Error![]HexCoords {
        var output: std.ArrayList(HexCoords) = .empty;
        defer output.deinit(allocator);
        for (self.hexes.items) |hex| for (hex.neighbors()) |neighbor| {
            if (self.contains(neighbor) or containsHex(output.items, neighbor)) continue;
            output.append(allocator, neighbor) catch return error.OutOfMemory;
        };
        return output.toOwnedSlice(allocator) catch return error.OutOfMemory;
    }

    pub fn isEquivalentWithout(self: Polyomino, coordinates: HexCoords) bool {
        if (self.countNeighbors(coordinates) != 3) return false;
        const neighbors = coordinates.neighbors();
        for (0..neighbors.len) |index| {
            if (self.contains(neighbors[index]) and self.contains(neighbors[(index + 5) % 6]) and self.contains(neighbors[(index + 4) % 6])) return true;
        }
        return false;
    }

    pub fn buildSkewedBox(self: *Polyomino, width: usize, height: usize, pentagon: bool) core.errors.Error!void {
        self.clear();
        for (0..height) |y| for (0..width) |x| try self.addHex(.{ .x = @intCast(x), .y = @intCast(y) });
        if (pentagon) try self.markOneVertexAsPentagon();
    }

    pub fn buildRaggedBox(self: *Polyomino, width: usize, height: usize, pentagon: bool) core.errors.Error!void {
        try self.buildRagged(width, height, .same, pentagon);
    }

    pub fn buildRaggedSmallerBox(self: *Polyomino, width: usize, height: usize, pentagon: bool) core.errors.Error!void {
        try self.buildRagged(width, height, .smaller, pentagon);
    }

    pub fn buildRaggedBiggerBox(self: *Polyomino, width: usize, height: usize, pentagon: bool) core.errors.Error!void {
        try self.buildRagged(width, height, .bigger, pentagon);
    }

    const RaggedKind = enum { same, smaller, bigger };

    fn buildRagged(self: *Polyomino, width: usize, height: usize, kind: RaggedKind, pentagon: bool) core.errors.Error!void {
        if (width == 0 or (kind == .smaller and width < 2)) return error.InvalidOption;
        self.clear();
        var start_x: i32 = 0;
        var y: usize = 0;
        while (y < height) {
            for (0..width) |x| try self.addHex(.{ .x = start_x + @as(i32, @intCast(x)), .y = @intCast(y) });
            y += 1;
            if (y >= height) break;
            const second_width = switch (kind) {
                .same => width,
                .smaller => width - 1,
                .bigger => width + 1,
            };
            const offset: i32 = if (kind == .bigger) -1 else 0;
            for (0..second_width) |x| try self.addHex(.{
                .x = start_x + offset + @as(i32, @intCast(x)),
                .y = @intCast(y),
            });
            start_x -= 1;
            y += 1;
        }
        if (pentagon) try self.markOneVertexAsPentagon();
    }

    pub fn buildWithVertices(self: *Polyomino, total_vertices: usize) core.errors.Error!void {
        if (total_vertices < 10) return error.InvalidOption;
        self.clear();
        try self.addHex(.{ .x = 0, .y = 0 });
        try self.addHex(.{ .x = 1, .y = 0 });
        var vertices: usize = 10;
        while (vertices < total_vertices) : (vertices += 2) {
            const candidates = try self.allFreeNeighbors(self.allocator);
            defer self.allocator.free(candidates);
            var best: ?usize = null;
            var lowest_distance: u32 = 0;
            for (candidates, 0..) |candidate, index| {
                if (self.countNeighbors(candidate) != 2) continue;
                const distance = candidate.distanceFrom(.{ .x = 0, .y = 0 });
                if (best == null or distance < lowest_distance) {
                    best = index;
                    lowest_distance = distance;
                }
            }
            const selected = best orelse return error.InvalidMapping;
            try self.addHex(candidates[selected]);
            for (candidates, 0..) |candidate, index| {
                if (index != selected and self.countNeighbors(candidate) == 3 and !self.contains(candidate)) try self.addHex(candidate);
            }
        }
        if (vertices - total_vertices == 1) try self.markOneVertexAsPentagon();
    }

    pub fn markOneVertexAsPentagon(self: *Polyomino) core.errors.Error!void {
        const contour = try self.path(self.allocator);
        defer self.allocator.free(contour);
        if (contour.len == 0) return error.InvalidMapping;
        const patterns = [_][3]usize{ .{ 2, 1, 2 }, .{ 1, 2, 1 } };
        for (patterns) |pattern| for (contour, 0..) |vertex, index| {
            const previous = try self.hexagonsAtVertex(contour[(index + contour.len - 1) % contour.len]);
            const current = try self.hexagonsAtVertex(vertex);
            const next = try self.hexagonsAtVertex(contour[(index + 1) % contour.len]);
            if (previous == pattern[0] and current == pattern[1] and next == pattern[2]) {
                self.pentagon_vertices.append(self.allocator, vertex) catch return error.OutOfMemory;
                return;
            }
        };
    }
};

pub const ShapeCollection = struct {
    allocator: std.mem.Allocator,
    items: std.ArrayList(Polyomino) = .empty,

    pub fn init(allocator: std.mem.Allocator) ShapeCollection {
        return .{ .allocator = allocator };
    }

    pub fn deinit(self: *ShapeCollection) void {
        for (self.items.items) |*item| item.deinit();
        self.items.deinit(self.allocator);
        self.* = undefined;
    }

    fn prepend(self: *ShapeCollection, shape: Polyomino) core.errors.Error!void {
        self.items.insert(self.allocator, 0, shape) catch return error.OutOfMemory;
    }

    pub fn removeDuplicates(self: *ShapeCollection) core.errors.Error!void {
        var index: usize = 0;
        while (index < self.items.items.len) {
            var duplicate = false;
            for (self.items.items[0..index]) |prior| {
                if (try self.items.items[index].isTheSameAs(prior, self.allocator)) {
                    duplicate = true;
                    break;
                }
            }
            if (duplicate) {
                var removed = self.items.orderedRemove(index);
                removed.deinit();
            } else index += 1;
        }
    }

    pub fn equivalents(self: ShapeCollection) core.errors.Error!ShapeCollection {
        var output = ShapeCollection.init(self.allocator);
        errdefer output.deinit();
        for (self.items.items) |shape| {
            const pentagons = shape.pentagon_vertices.items.len;
            for (shape.hexes.items) |coordinates| {
                if (!shape.isEquivalentWithout(coordinates)) continue;
                var equivalent = try shape.clone(self.allocator);
                errdefer equivalent.deinit();
                equivalent.pentagon_vertices.clearRetainingCapacity();
                try equivalent.removeHex(coordinates);
                for (0..pentagons) |_| try equivalent.markOneVertexAsPentagon();
                output.items.append(self.allocator, equivalent) catch return error.OutOfMemory;
            }
        }
        return output;
    }
};

const ShapeKind = enum { skewed, ragged, ragged_smaller, ragged_bigger };

fn prependShape(
    collection: *ShapeCollection,
    kind: ShapeKind,
    width: usize,
    height: usize,
    pentagon: bool,
) core.errors.Error!void {
    var shape = Polyomino.init(collection.allocator);
    errdefer shape.deinit();
    switch (kind) {
        .skewed => try shape.buildSkewedBox(width, height, pentagon),
        .ragged => try shape.buildRaggedBox(width, height, pentagon),
        .ragged_smaller => try shape.buildRaggedSmallerBox(width, height, pentagon),
        .ragged_bigger => try shape.buildRaggedBiggerBox(width, height, pentagon),
    }
    try collection.prepend(shape);
}

pub fn buildSquaredShapes(allocator: std.mem.Allocator, requested_vertices: usize) core.errors.Error!ShapeCollection {
    if (requested_vertices < 3 or requested_vertices == std.math.maxInt(usize)) return error.InvalidOption;
    var output = ShapeCollection.init(allocator);
    errdefer output.deinit();
    const pentagon = requested_vertices % 2 != 0;
    const total_vertices = requested_vertices + @intFromBool(pentagon);
    if (total_vertices % 4 == 0) {
        if (total_vertices >= 12) {
            const sum = total_vertices / 4;
            const middle = sum / 2;
            for (1..middle) |first| {
                const second = sum - first;
                if (second % 2 == 0 and first > 1) try prependShape(&output, .ragged_bigger, first, second, pentagon);
                if (first % 2 == 0 and second > 1) try prependShape(&output, .ragged_bigger, second, first, pentagon);
            }
        }
    } else {
        const sum = (total_vertices + 2) / 4;
        const middle = sum / 2;
        for (1..middle + 1) |first| {
            const second = sum - first;
            try prependShape(&output, .skewed, first, second, pentagon);
            if (first < 2 or second < 2) continue;
            try prependShape(&output, .ragged, first, second, pentagon);
            try prependShape(&output, .ragged, second, first, pentagon);
            if (second % 2 != 0) try prependShape(&output, .ragged_bigger, first, second, pentagon);
            if (first % 2 != 0) try prependShape(&output, .ragged_bigger, second, first, pentagon);
            if (first > 2 and second % 2 != 0) try prependShape(&output, .ragged_smaller, first, second, pentagon);
            if (second > 2 and first % 2 != 0) try prependShape(&output, .ragged_smaller, second, first, pentagon);
        }
    }
    return output;
}

pub const DoubleBondConstraint = struct {
    trans: bool,
    previous_atom: usize,
    atom_a: usize,
    atom_b: usize,
    following_atom: usize,
};

pub const RingConstraint = struct {
    ring: core.ids.RingId,
    atom: usize,
    force_outside: bool,
};

pub const SubstitutedAtom = struct {
    atom: usize,
    subtree_size: usize,
};

pub const PathConstraints = struct {
    double_bonds: []const DoubleBondConstraint = &.{},
    rings: []const RingConstraint = &.{},
};

pub const PathRestraints = struct {
    hetero_atoms: []const usize = &.{},
    substituted_atoms: []const SubstitutedAtom = &.{},
};

pub const Match = struct {
    start: usize,
    score: i32,
};

pub const ShapeMatch = struct {
    shape: usize,
    start: usize,
    score: i32,
};

pub fn lowestPeriod(neighbors: []const u8) usize {
    for (1..neighbors.len) |period| {
        var different = false;
        for (neighbors, 0..) |value, index| {
            if (value != neighbors[(index + period) % neighbors.len]) {
                different = true;
                break;
            }
        }
        if (!different) return period;
    }
    return neighbors.len;
}

pub fn matchPolyomino(
    allocator: std.mem.Allocator,
    shape: Polyomino,
    constraints: PathConstraints,
    restraints: PathRestraints,
) core.errors.Error!?Match {
    const contour = try shape.path(allocator);
    defer allocator.free(contour);
    if (contour.len == 0) return null;
    const neighbors = allocator.alloc(u8, contour.len) catch return error.OutOfMemory;
    defer allocator.free(neighbors);
    for (contour, neighbors) |vertex, *count| count.* = @intCast(try shape.hexagonsAtVertex(vertex));
    var best: ?Match = null;
    for (0..lowestPeriod(neighbors)) |start| {
        const score = try scorePath(allocator, shape, contour, neighbors, start, constraints, restraints) orelse continue;
        if (best == null or score > best.?.score) {
            best = .{ .start = start, .score = score };
            if (score == 0) break;
        }
    }
    return best;
}

pub fn matchShapes(
    allocator: std.mem.Allocator,
    shapes: ShapeCollection,
    constraints: PathConstraints,
    restraints: PathRestraints,
    checked_macrocycles: *usize,
) core.errors.Error!?ShapeMatch {
    var best: ?ShapeMatch = null;
    for (shapes.items.items, 0..) |shape, index| {
        if (try matchPolyomino(allocator, shape, constraints, restraints)) |matched| {
            if (best == null or matched.score > best.?.score) {
                best = .{ .shape = index, .start = matched.start, .score = matched.score };
                if (matched.score == 0) return best;
            }
        }
        const previous = checked_macrocycles.*;
        checked_macrocycles.* += 1;
        if (previous > max_macrocycles) break;
    }
    return best;
}

fn scorePath(
    allocator: std.mem.Allocator,
    shape: Polyomino,
    contour: []const VertexCoords,
    neighbors: []const u8,
    start: usize,
    constraints: PathConstraints,
    restraints: PathRestraints,
) core.errors.Error!?i32 {
    if (!try checkRingConstraints(allocator, shape, contour, neighbors, start, constraints.rings)) return null;
    if (!checkDoubleBondConstraints(contour, start, constraints.double_bonds)) return null;
    var score: i32 = 0;
    for (restraints.hetero_atoms) |atom| {
        if (atom >= contour.len) return error.InvalidMapping;
        if (neighbors[(atom + start) % contour.len] == 1) score -= 1;
    }
    var used_substituents: std.ArrayList(VertexCoords) = .empty;
    defer used_substituents.deinit(allocator);
    for (restraints.substituted_atoms, 0..) |substituted, index| {
        if (substituted.atom >= contour.len or index >= contour.len or substituted.subtree_size > @as(usize, std.math.maxInt(i32) / 10)) return error.InvalidMapping;
        const counter = (substituted.atom + start) % contour.len;
        if (neighbors[counter] != 2) continue;
        score -= @as(i32, @intCast(substituted.subtree_size)) * 10;
        // Preserve upstream's observable indexing bug: restraint index rather
        // than mapped atom index selects the substituent projection vertex.
        const projected = try shape.coordinatesOfSubstituent(contour[index]);
        if (containsVertex(used_substituents.items, projected)) score -= 200;
        if (containsVertex(contour, projected)) score -= 400;
        used_substituents.append(allocator, projected) catch return error.OutOfMemory;
    }
    return score;
}

fn checkDoubleBondConstraints(
    contour: []const VertexCoords,
    start: usize,
    constraints: []const DoubleBondConstraint,
) bool {
    for (constraints) |constraint| {
        if (constraint.previous_atom >= contour.len or constraint.atom_a >= contour.len or
            constraint.atom_b >= contour.len or constraint.following_atom >= contour.len) return false;
        const previous = contour[(start + constraint.previous_atom) % contour.len].toCartesian();
        const atom_a = contour[(start + constraint.atom_a) % contour.len].toCartesian();
        const atom_b = contour[(start + constraint.atom_b) % contour.len].toCartesian();
        const following = contour[(start + constraint.following_atom) % contour.len].toCartesian();
        if (geometry.sameSide(previous, following, atom_a, atom_b) == constraint.trans) return false;
    }
    return true;
}

const AllowedRing = struct {
    ring: core.ids.RingId,
    positions: [3]HexCoords,
    len: usize,
};

fn checkRingConstraints(
    allocator: std.mem.Allocator,
    shape: Polyomino,
    contour: []const VertexCoords,
    neighbors: []const u8,
    start: usize,
    constraints: []const RingConstraint,
) core.errors.Error!bool {
    var allowed: std.ArrayList(AllowedRing) = .empty;
    defer allowed.deinit(allocator);
    for (constraints) |constraint| {
        if (constraint.atom >= contour.len or !constraint.ring.isValid()) return error.InvalidMapping;
        const counter = (constraint.atom + start) % contour.len;
        if (constraint.force_outside and neighbors[counter] != 1) return false;
        var positions: [3]HexCoords = undefined;
        const position_count = try shape.freeVertexNeighborPositions(contour[counter], &positions);
        var existing_index: ?usize = null;
        for (allowed.items, 0..) |entry, index| if (entry.ring == constraint.ring) {
            existing_index = index;
            break;
        };
        if (existing_index) |index| {
            var intersection: [3]HexCoords = undefined;
            var count: usize = 0;
            for (positions[0..position_count]) |position| {
                if (containsHex(allowed.items[index].positions[0..allowed.items[index].len], position)) {
                    intersection[count] = position;
                    count += 1;
                }
            }
            if (count == 0) return false;
            allowed.items[index].positions = intersection;
            allowed.items[index].len = count;
        } else {
            if (position_count == 0) return false;
            allowed.append(allocator, .{ .ring = constraint.ring, .positions = positions, .len = position_count }) catch return error.OutOfMemory;
        }
    }
    return true;
}

pub fn writeCoordinates(
    shape: Polyomino,
    start: usize,
    coordinates: []core.math.Vec2,
    coordinates_set: []const bool,
) core.errors.Error!void {
    if (coordinates.len != coordinates_set.len) return error.InvalidMapping;
    const contour = try shape.path(shape.allocator);
    defer shape.allocator.free(contour);
    if (coordinates.len != contour.len or start >= contour.len) return error.InvalidMapping;
    for (coordinates, coordinates_set, 0..) |*coordinate, already_set, index| {
        if (!already_set) coordinate.* = contour[(index + start) % contour.len].toCartesian();
    }
}

fn containsHex(items: []const HexCoords, target: HexCoords) bool {
    for (items) |item| if (std.meta.eql(item, target)) return true;
    return false;
}

fn containsVertex(items: []const VertexCoords, target: VertexCoords) bool {
    for (items) |item| if (std.meta.eql(item, target)) return true;
    return false;
}

fn makePolyomino(allocator: std.mem.Allocator, coordinates: []const HexCoords) !Polyomino {
    var result = Polyomino.init(allocator);
    errdefer result.deinit();
    for (coordinates) |coordinate| try result.addHex(coordinate);
    return result;
}

test "polyomino substituent coordinates match rehosted upstream cases" {
    var polyomino = Polyomino.init(std.testing.allocator);
    defer polyomino.deinit();
    try polyomino.addHex(.{ .x = 0, .y = 0 });
    try std.testing.expectEqual(VertexCoords{ .x = 1, .y = -1, .z = -1 }, try polyomino.coordinatesOfSubstituent(.{ .x = 1, .y = 0, .z = 0 }));
    try polyomino.addHex(.{ .x = 1, .y = 0 });
    try std.testing.expectEqual(VertexCoords{ .x = 0, .y = 0, .z = -1 }, try polyomino.coordinatesOfSubstituent(.{ .x = 1, .y = 0, .z = 0 }));
}

test "polyomino identity supports rotations but rejects reflection" {
    var first = try makePolyomino(std.testing.allocator, &.{ .{ .x = 0, .y = 0 }, .{ .x = 1, .y = 0 }, .{ .x = 2, .y = 0 }, .{ .x = 0, .y = 1 } });
    defer first.deinit();
    var translated = try makePolyomino(std.testing.allocator, &.{ .{ .x = 4, .y = 2 }, .{ .x = 5, .y = 2 }, .{ .x = 6, .y = 2 }, .{ .x = 4, .y = 3 } });
    defer translated.deinit();
    var rotated = try makePolyomino(std.testing.allocator, &.{ .{ .x = 0, .y = 0 }, .{ .x = -1, .y = 0 }, .{ .x = -2, .y = 0 }, .{ .x = 0, .y = -1 } });
    defer rotated.deinit();
    var reflected = try makePolyomino(std.testing.allocator, &.{ .{ .x = 0, .y = 0 }, .{ .x = 0, .y = 1 }, .{ .x = 0, .y = 2 }, .{ .x = 1, .y = 0 } });
    defer reflected.deinit();
    try std.testing.expect(try first.isTheSameAs(translated, std.testing.allocator));
    try std.testing.expect(try translated.isTheSameAs(first, std.testing.allocator));
    try std.testing.expect(try first.isTheSameAs(rotated, std.testing.allocator));
    try std.testing.expect(!try first.isTheSameAs(reflected, std.testing.allocator));
}

test "round polyomino path is deterministic and allocation failures clean up" {
    var polyomino = Polyomino.init(std.testing.allocator);
    defer polyomino.deinit();
    try polyomino.buildWithVertices(17);
    const contour = try polyomino.path(std.testing.allocator);
    defer std.testing.allocator.free(contour);
    try std.testing.expectEqual(@as(usize, 17), contour.len);
}

fn buildAndTraverse(allocator: std.mem.Allocator) !void {
    var polyomino = Polyomino.init(allocator);
    defer polyomino.deinit();
    try polyomino.buildWithVertices(17);
    const contour = try polyomino.path(allocator);
    defer allocator.free(contour);
    var copy = try polyomino.clone(allocator);
    defer copy.deinit();
    try std.testing.expect(try polyomino.isTheSameAs(copy, allocator));
}

test "polyomino allocation failures leave no owned lattice state" {
    try std.testing.checkAllAllocationFailures(std.testing.allocator, buildAndTraverse, .{});
}

test "squared shape enumeration matches pinned candidate and path order" {
    const cases = [_]struct { vertices: usize, hexes: []const usize }{
        .{ .vertices = 9, .hexes = &.{2} },
        .{ .vertices = 11, .hexes = &.{} },
        .{ .vertices = 13, .hexes = &.{ 4, 4, 4, 3 } },
        .{ .vertices = 17, .hexes = &.{ 7, 6, 6, 6, 4 } },
        .{ .vertices = 21, .hexes = &.{ 8, 8, 10, 10, 9, 9, 9, 8, 8, 8, 5 } },
    };
    for (cases) |case| {
        var shapes = try buildSquaredShapes(std.testing.allocator, case.vertices);
        defer shapes.deinit();
        try std.testing.expectEqual(case.hexes.len, shapes.items.items.len);
        for (shapes.items.items, case.hexes) |shape, expected_hexes| {
            try std.testing.expectEqual(expected_hexes, shape.hexes.items.len);
            const contour = try shape.path(std.testing.allocator);
            defer std.testing.allocator.free(contour);
            try std.testing.expectEqual(case.vertices, contour.len);
        }
    }
}

fn enumerateAndDeduplicate(allocator: std.mem.Allocator) !void {
    var shapes = try buildSquaredShapes(allocator, 21);
    defer shapes.deinit();
    try shapes.removeDuplicates();
    var equivalents = try shapes.equivalents();
    defer equivalents.deinit();
    try equivalents.removeDuplicates();
}

test "shape enumeration equivalence and deduplication clean every allocation failure" {
    try std.testing.checkAllAllocationFailures(std.testing.allocator, enumerateAndDeduplicate, .{});
}

fn matchAndWrite(allocator: std.mem.Allocator) !void {
    var shape = Polyomino.init(allocator);
    defer shape.deinit();
    try shape.buildWithVertices(17);
    const contour = try shape.path(allocator);
    defer allocator.free(contour);
    var inward: ?usize = null;
    for (contour, 0..) |vertex, index| {
        if (try shape.hexagonsAtVertex(vertex) == 2) {
            inward = index;
            break;
        }
    }
    const substituted = [_]SubstitutedAtom{.{ .atom = inward orelse return error.InvalidMapping, .subtree_size = 2 }};
    const rings = [_]RingConstraint{.{ .ring = core.ids.RingId.fromIndex(0), .atom = 0, .force_outside = false }};
    const matched = (try matchPolyomino(allocator, shape, .{ .rings = &rings }, .{ .substituted_atoms = &substituted })) orelse return error.InvalidMapping;
    var coordinates = try allocator.alloc(core.math.Vec2, contour.len);
    defer allocator.free(coordinates);
    @memset(coordinates, .{});
    const coordinates_set = try allocator.alloc(bool, contour.len);
    defer allocator.free(coordinates_set);
    @memset(coordinates_set, false);
    coordinates_set[0] = true;
    coordinates[0] = .{ .x = 123, .y = 456 };
    try writeCoordinates(shape, matched.start, coordinates, coordinates_set);
    try std.testing.expectEqual(core.math.Vec2{ .x = 123, .y = 456 }, coordinates[0]);
    try std.testing.expect(coordinates[1].isFinite());
}

test "path matching writes only unset coordinates and cleans allocation failures" {
    try std.testing.checkAllAllocationFailures(std.testing.allocator, matchAndWrite, .{});
}

test "invalid hard constraints reject every rotational start" {
    var shape = Polyomino.init(std.testing.allocator);
    defer shape.deinit();
    try shape.buildWithVertices(17);
    const invalid = [_]DoubleBondConstraint{.{
        .trans = true,
        .previous_atom = 99,
        .atom_a = 0,
        .atom_b = 1,
        .following_atom = 2,
    }};
    try std.testing.expect((try matchPolyomino(std.testing.allocator, shape, .{ .double_bonds = &invalid }, .{})) == null);
}
