const std = @import("std");
const core = @import("core");

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

    pub fn buildSkewedBox(self: *Polyomino, width: usize, height: usize) core.errors.Error!void {
        self.clear();
        for (0..height) |y| for (0..width) |x| try self.addHex(.{ .x = @intCast(x), .y = @intCast(y) });
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
