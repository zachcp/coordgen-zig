const std = @import("std");
const core = @import("core");

pub const Contour = struct {
    allocator: std.mem.Allocator,
    points: []core.math.Vec2,
    edges: []Edge,
    visited: []bool,

    pub const Edge = struct { start: u32, end: u32 };

    pub fn deinit(self: *Contour) void {
        self.allocator.free(self.visited);
        self.allocator.free(self.edges);
        self.allocator.free(self.points);
        self.* = undefined;
    }

    /// Return independent coordinate copies in upstream's destructive
    /// point-allocation traversal order. A second call returns no contours.
    pub fn orderedCoordinates(self: *Contour, allocator: std.mem.Allocator) core.errors.Error!OrderedContours {
        var points: std.ArrayList(core.math.Vec2) = .empty;
        defer points.deinit(allocator);
        var offsets: std.ArrayList(u32) = .empty;
        defer offsets.deinit(allocator);
        offsets.append(allocator, 0) catch return error.OutOfMemory;

        while (std.mem.indexOfScalar(bool, self.visited, false)) |start| {
            var current: ?u32 = @intCast(start);
            while (current) |index| {
                self.visited[index] = true;
                points.append(allocator, self.points[index]) catch return error.OutOfMemory;
                current = self.firstUnvisitedNeighbor(index);
            }
            if (points.items.len > std.math.maxInt(u32)) return error.TooManyItems;
            offsets.append(allocator, @intCast(points.items.len)) catch return error.OutOfMemory;
        }

        const owned_points = points.toOwnedSlice(allocator) catch return error.OutOfMemory;
        errdefer allocator.free(owned_points);
        const owned_offsets = offsets.toOwnedSlice(allocator) catch return error.OutOfMemory;
        return .{ .allocator = allocator, .points = owned_points, .offsets = owned_offsets };
    }

    fn firstUnvisitedNeighbor(self: Contour, point: u32) ?u32 {
        var found: usize = 0;
        for (self.edges) |edge| {
            const neighbor = if (edge.start == point)
                edge.end
            else if (edge.end == point)
                edge.start
            else
                continue;
            found += 1;
            if (!self.visited[neighbor]) return neighbor;
            if (found == 2) return null;
        }
        return null;
    }
};

pub const OrderedContours = struct {
    allocator: std.mem.Allocator,
    points: []core.math.Vec2,
    offsets: []u32,

    pub fn deinit(self: *OrderedContours) void {
        self.allocator.free(self.offsets);
        self.allocator.free(self.points);
        self.* = undefined;
    }

    pub fn count(self: OrderedContours) usize {
        return self.offsets.len - 1;
    }

    pub fn contour(self: OrderedContours, index: usize) []const core.math.Vec2 {
        return self.points[self.offsets[index]..self.offsets[index + 1]];
    }
};

/// Deterministic bottom-to-top, left-to-right scalar contour extraction. Edge
/// reuse and checkerboard disambiguation follow the pinned upstream sweep.
pub fn run(
    allocator: std.mem.Allocator,
    values: []const f32,
    width: usize,
    height: usize,
    origin: core.math.Vec2,
    interval: core.math.Vec2,
    threshold: f32,
) core.errors.Error!Contour {
    if (width < 2 or height < 2 or values.len != width * height) return error.InvalidMapping;
    if (!origin.isFinite() or !interval.isFinite() or !std.math.isFinite(threshold) or interval.x <= 0 or interval.y <= 0) return error.InvalidOption;
    for (values) |value| if (!std.math.isFinite(value)) return error.InvalidCoordinate;
    var points: std.ArrayList(core.math.Vec2) = .empty;
    defer points.deinit(allocator);
    var edges: std.ArrayList(Contour.Edge) = .empty;
    defer edges.deinit(allocator);
    const previous_top = allocator.alloc(?u32, width - 1) catch return error.OutOfMemory;
    defer allocator.free(previous_top);
    @memset(previous_top, null);

    for (0..height - 1) |y| {
        var previous_right: ?u32 = null;
        for (0..width - 1) |x| {
            const bottom_left = values[x + y * width];
            const bottom_right = values[x + 1 + y * width];
            const top_left = values[x + (y + 1) * width];
            const top_right = values[x + 1 + (y + 1) * width];
            const left = previous_right;
            const bottom = previous_top[x];
            const right = try crossing(allocator, &points, .{ .x = @floatFromInt(x + 1), .y = @floatFromInt(y) }, .{ .x = @floatFromInt(x + 1), .y = @floatFromInt(y + 1) }, bottom_right, top_right, origin, interval, threshold);
            const top = try crossing(allocator, &points, .{ .x = @floatFromInt(x), .y = @floatFromInt(y + 1) }, .{ .x = @floatFromInt(x + 1), .y = @floatFromInt(y + 1) }, top_left, top_right, origin, interval, threshold);
            previous_right = right;
            previous_top[x] = top;
            if (top != null and right != null and bottom != null and left != null) {
                if (top_left > threshold) {
                    try appendEdge(allocator, &edges, top.?, right.?);
                    try appendEdge(allocator, &edges, bottom.?, left.?);
                } else {
                    try appendEdge(allocator, &edges, top.?, left.?);
                    try appendEdge(allocator, &edges, bottom.?, right.?);
                }
                continue;
            }
            const candidates = [_]?u32{ top, right, bottom, left };
            var selected: [2]u32 = undefined;
            var count: usize = 0;
            for (candidates) |candidate| if (candidate) |point| {
                if (count == 0) {
                    selected[0] = point;
                } else {
                    // Upstream repeatedly overwrites p2, so malformed or
                    // boundary three-crossing cells connect first to last.
                    selected[1] = point;
                }
                count += 1;
            };
            if (count >= 2) try appendEdge(allocator, &edges, selected[0], selected[1]);
            _ = bottom_left;
        }
    }
    const owned_points = points.toOwnedSlice(allocator) catch return error.OutOfMemory;
    errdefer allocator.free(owned_points);
    const owned_edges = edges.toOwnedSlice(allocator) catch return error.OutOfMemory;
    errdefer allocator.free(owned_edges);
    const visited = allocator.alloc(bool, owned_points.len) catch return error.OutOfMemory;
    @memset(visited, false);
    return .{ .allocator = allocator, .points = owned_points, .edges = owned_edges, .visited = visited };
}

fn crossing(
    allocator: std.mem.Allocator,
    points: *std.ArrayList(core.math.Vec2),
    first: core.math.Vec2,
    second: core.math.Vec2,
    first_value: f32,
    second_value: f32,
    origin: core.math.Vec2,
    interval: core.math.Vec2,
    threshold: f32,
) core.errors.Error!?u32 {
    if ((first_value - threshold) * (second_value - threshold) >= 0) return null;
    const difference = second_value - first_value;
    const parameter: f32 = if (@abs(difference) < 0.0001) 0.5 else (threshold - first_value) / difference;
    const grid = core.math.Vec2{
        .x = first.x + (second.x - first.x) * parameter,
        .y = first.y + (second.y - first.y) * parameter,
    };
    if (points.items.len >= std.math.maxInt(u32)) return error.TooManyItems;
    const index: u32 = @intCast(points.items.len);
    points.append(allocator, .{
        .x = origin.x + grid.x * interval.x,
        .y = origin.y + grid.y * interval.y,
    }) catch return error.OutOfMemory;
    return index;
}

fn appendEdge(allocator: std.mem.Allocator, edges: *std.ArrayList(Contour.Edge), start: u32, end: u32) core.errors.Error!void {
    edges.append(allocator, .{ .start = start, .end = end }) catch return error.OutOfMemory;
}

fn runAndDiscard(allocator: std.mem.Allocator, values: []const f32) !void {
    var contour = try run(allocator, values, 3, 3, .{}, .{ .x = 1, .y = 1 }, 0.5);
    contour.deinit();
}

fn runOrderAndDiscard(allocator: std.mem.Allocator, values: []const f32) !void {
    var contour = try run(allocator, values, 3, 3, .{}, .{ .x = 1, .y = 1 }, 0.5);
    defer contour.deinit();
    var ordered = try contour.orderedCoordinates(allocator);
    ordered.deinit();
}

test "marching squares has deterministic point creation and checkerboard topology" {
    const values = [_]f32{ 0, 1, 0, 1, 0, 1, 0, 1, 0 };
    var first = try run(std.testing.allocator, &values, 3, 3, .{}, .{ .x = 10, .y = 20 }, 0.5);
    defer first.deinit();
    var second = try run(std.testing.allocator, &values, 3, 3, .{}, .{ .x = 10, .y = 20 }, 0.5);
    defer second.deinit();
    try std.testing.expectEqualSlices(core.math.Vec2, first.points, second.points);
    try std.testing.expectEqualSlices(Contour.Edge, first.edges, second.edges);
    try std.testing.expectEqualSlices(Contour.Edge, &.{
        .{ .start = 1, .end = 0 },
        .{ .start = 3, .end = 0 },
        .{ .start = 5, .end = 1 },
        .{ .start = 7, .end = 6 },
        .{ .start = 3, .end = 4 },
    }, first.edges);

    var ordered = try first.orderedCoordinates(std.testing.allocator);
    defer ordered.deinit();
    try std.testing.expectEqual(@as(usize, 4), ordered.count());
    try std.testing.expectEqualSlices(core.math.Vec2, &.{
        .{ .x = 10, .y = 10 },
        .{ .x = 5, .y = 20 },
        .{ .x = 5, .y = 40 },
    }, ordered.contour(0));
    try std.testing.expectEqualSlices(core.math.Vec2, &.{.{ .x = 20, .y = 10 }}, ordered.contour(1));
    var exhausted = try first.orderedCoordinates(std.testing.allocator);
    defer exhausted.deinit();
    try std.testing.expectEqual(@as(usize, 0), exhausted.count());
}

test "marching squares rejects nonfinite grids and cleans allocation failures" {
    try std.testing.expectError(error.InvalidCoordinate, run(std.testing.allocator, &.{ 0, 1, 0, std.math.nan(f32) }, 2, 2, .{}, .{ .x = 1, .y = 1 }, 0.5));
    const values = [_]f32{ 0, 1, 0, 1, 0, 1, 0, 1, 0 };
    try std.testing.checkAllAllocationFailures(std.testing.allocator, runAndDiscard, .{&values});
    try std.testing.checkAllAllocationFailures(std.testing.allocator, runOrderAndDiscard, .{&values});
}
