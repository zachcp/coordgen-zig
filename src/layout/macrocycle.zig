const std = @import("std");
const core = @import("core");
const geometry = @import("geometry");
const model = @import("model");
const topology = @import("topology");
const fragments = @import("fragments.zig");
const neighbour_order = @import("neighbour_order.zig");

pub const max_macrocycles: usize = 40;
pub const path_failed: i32 = -1000;
pub const substituted_atom_restraint: i32 = 10;
pub const sqrt_three_halves: f32 = 0.8660254037844386;

/// Pinned `CoordgenMacrocycleBuilder::acceptableShapeScore`. Path scores start
/// at zero and only ever subtract, so for ten atoms or more the threshold is
/// unreachable and the early break never fires: upstream keeps expanding the
/// shape pool until the `max_macrocycles` cap or a perfect zero score. The
/// positive sign is upstream's and is load-bearing, not a transcription slip.
pub fn acceptableShapeScore(atom_count: usize) core.errors.Error!i32 {
    if (atom_count < 10) return 0;
    if (atom_count > @as(usize, @intCast(@divTrunc(std.math.maxInt(i32), substituted_atom_restraint)))) return error.TooManyItems;
    return @divTrunc(@as(i32, @intCast(atom_count)) * substituted_atom_restraint, 2);
}

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
        // Two adjacent hexagons start with ten contour vertices; replacing
        // one shared-corner vertex with a pentagon yields the valid nine-atom
        // macrocycle shape used at upstream's minimum macrocycle size.
        if (total_vertices < 9) return error.InvalidOption;
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

pub const PathData = struct {
    allocator: std.mem.Allocator,
    ordered_atoms: []core.ids.AtomId,
    double_bonds: []DoubleBondConstraint,
    rings: []RingConstraint,
    hetero_atoms: []usize,
    substituted_atoms: []SubstitutedAtom,

    pub fn deinit(self: *PathData) void {
        self.allocator.free(self.substituted_atoms);
        self.allocator.free(self.hetero_atoms);
        self.allocator.free(self.rings);
        self.allocator.free(self.double_bonds);
        self.allocator.free(self.ordered_atoms);
        self.* = undefined;
    }

    pub fn constraints(self: PathData) PathConstraints {
        return .{ .double_bonds = self.double_bonds, .rings = self.rings };
    }

    pub fn restraints(self: PathData) PathRestraints {
        return .{ .hetero_atoms = self.hetero_atoms, .substituted_atoms = self.substituted_atoms };
    }
};

pub const GenerateResult = enum {
    matched,
    /// A shape was applied, and it carries pentagon vertices. Upstream calls
    /// requireMinimization in exactly this case
    /// (CoordgenMacrocycleBuilder.cpp, after writePolyominoCoordinates),
    /// because a pentagon vertex leaves the ring off the hexagonal lattice
    /// and the geometry needs relaxing.
    matched_needs_minimization,
    no_shape,

    pub fn placed(self: GenerateResult) bool {
        return self != .no_shape;
    }
};

/// Extract the ordered path and every hard/soft macrocycle condition from the
/// native topology. The owned result keeps no pointers into temporary search
/// state and can be reused across equivalent-shape rounds.
pub fn collectPathData(
    allocator: std.mem.Allocator,
    ring: core.ids.RingId,
    atoms: []const model.Atom,
    bonds: []const model.Bond,
    graph: topology.Graph,
    membership: topology.RingMembership,
) core.errors.Error!PathData {
    const ordered = try orderRing(allocator, ring, graph, membership);
    errdefer allocator.free(ordered);
    const positions = allocator.alloc(usize, atoms.len) catch return error.OutOfMemory;
    defer allocator.free(positions);
    @memset(positions, std.math.maxInt(usize));
    for (ordered, 0..) |atom, index| positions[atom.index()] = index;

    var double_bonds: std.ArrayList(DoubleBondConstraint) = .empty;
    defer double_bonds.deinit(allocator);
    var ring_constraints: std.ArrayList(RingConstraint) = .empty;
    defer ring_constraints.deinit(allocator);
    var hetero_atoms: std.ArrayList(usize) = .empty;
    defer hetero_atoms.deinit(allocator);
    var substituted_atoms: std.ArrayList(SubstitutedAtom) = .empty;
    defer substituted_atoms.deinit(allocator);

    for (ordered, 0..) |atom_id, index| {
        if (atoms[atom_id.index()].atomic_number != .carbon) {
            hetero_atoms.append(allocator, index) catch return error.OutOfMemory;
        }
        if (graph.degree(atom_id) != 2) {
            var subtree_size: usize = 0;
            const previous = ordered[(index + ordered.len - 1) % ordered.len];
            const following = ordered[(index + 1) % ordered.len];
            for (graph.neighbors(atom_id)) |neighbor| {
                if (neighbor == previous or neighbor == following) continue;
                const subtree = try graph.reachableExcluding(allocator, neighbor, atom_id);
                defer allocator.free(subtree);
                subtree_size = std.math.add(usize, subtree_size, subtree.len) catch return error.TooManyItems;
            }
            substituted_atoms.append(allocator, .{ .atom = index, .subtree_size = subtree_size }) catch return error.OutOfMemory;
        }

        if (membership.atomRings(atom_id).len > 1) {
            for (membership.atomRings(atom_id)) |other_ring| {
                if (membership.atoms(other_ring).len >= topology.rings.macrocycle_size) continue;
                var force_outside = false;
                for (graph.neighbors(atom_id)) |neighbor| {
                    if (positions[neighbor.index()] != std.math.maxInt(usize)) continue;
                    force_outside = std.mem.indexOfScalar(core.ids.AtomId, membership.atoms(other_ring), neighbor) != null;
                    break;
                }
                ring_constraints.append(allocator, .{
                    .ring = other_ring,
                    .atom = index,
                    .force_outside = force_outside,
                }) catch return error.OutOfMemory;
            }
        }

        const next_index = (index + 1) % ordered.len;
        const bond_id = bondBetween(atom_id, ordered[next_index], graph) orelse return error.InvalidMapping;
        const bond = bonds[bond_id.index()];
        if (bond.effective_order != .double) continue;
        var small_ring_bond = false;
        if (membership.bondRings(bond_id).len > 1) for (membership.bondRings(bond_id)) |bond_ring| {
            if (membership.atoms(bond_ring).len < topology.rings.macrocycle_size) {
                small_ring_bond = true;
                break;
            }
        };
        if (small_ring_bond) continue;
        const absolute = try topology.stereo.absoluteBondStereo(allocator, atoms, bonds, graph, membership, bond_id);
        if (absolute != .z and absolute != .e) continue;
        var previous_index = (index + ordered.len - 1) % ordered.len;
        var following_index = (index + 2) % ordered.len;
        var atom_a = index;
        var atom_b = next_index;
        if (bond.start != atom_id) {
            std.mem.swap(usize, &previous_index, &following_index);
            atom_a = next_index;
            atom_b = index;
        }
        var trans = absolute == .e;
        const start_first = (try topology.stereo.firstNeighbor(allocator, atoms, bonds, graph, bond.start, bond.end)) orelse continue;
        const end_first = (try topology.stereo.firstNeighbor(allocator, atoms, bonds, graph, bond.end, bond.start)) orelse continue;
        if (start_first != ordered[previous_index]) trans = !trans;
        if (end_first != ordered[following_index]) trans = !trans;
        double_bonds.append(allocator, .{
            .trans = trans,
            .previous_atom = previous_index,
            .atom_a = atom_a,
            .atom_b = atom_b,
            .following_atom = following_index,
        }) catch return error.OutOfMemory;
    }

    const owned_double_bonds = double_bonds.toOwnedSlice(allocator) catch return error.OutOfMemory;
    errdefer allocator.free(owned_double_bonds);
    const owned_rings = ring_constraints.toOwnedSlice(allocator) catch return error.OutOfMemory;
    errdefer allocator.free(owned_rings);
    const owned_hetero_atoms = hetero_atoms.toOwnedSlice(allocator) catch return error.OutOfMemory;
    errdefer allocator.free(owned_hetero_atoms);
    const owned_substituted_atoms = substituted_atoms.toOwnedSlice(allocator) catch return error.OutOfMemory;
    return .{
        .allocator = allocator,
        .ordered_atoms = ordered,
        .double_bonds = owned_double_bonds,
        .rings = owned_rings,
        .hetero_atoms = owned_hetero_atoms,
        .substituted_atoms = owned_substituted_atoms,
    };
}

/// Run the pinned round/squared/equivalent search and write the best valid
/// lattice path. `no_shape` deliberately leaves all coordinates untouched so
/// the caller can choose the open-cycle or regular-ring fallback.
pub fn generateShape(
    allocator: std.mem.Allocator,
    data: PathData,
    coordinates: []core.math.Vec2,
    coordinates_set: []const bool,
    force_open_macrocycles: bool,
) core.errors.Error!GenerateResult {
    if (coordinates.len != data.ordered_atoms.len or coordinates_set.len != coordinates.len) return error.InvalidMapping;
    var shapes = ShapeCollection.init(allocator);
    defer shapes.deinit();
    var round = Polyomino.init(allocator);
    errdefer round.deinit();
    try round.buildWithVertices(data.ordered_atoms.len);
    shapes.items.append(allocator, round) catch return error.OutOfMemory;
    round = Polyomino.init(allocator);
    var squared = try buildSquaredShapes(allocator, data.ordered_atoms.len);
    defer squared.deinit();
    while (squared.items.items.len != 0) {
        var shape = squared.items.orderedRemove(0);
        errdefer shape.deinit();
        shapes.items.append(allocator, shape) catch return error.OutOfMemory;
        shape = Polyomino.init(allocator);
    }
    try shapes.removeDuplicates();
    if (force_open_macrocycles) return .no_shape;

    var checked: usize = 0;
    var chosen: ?Polyomino = null;
    defer {
        if (chosen) |*shape| shape.deinit();
    }
    var chosen_start: usize = 0;
    var chosen_score: i32 = path_failed;
    var found = false;
    const acceptable_score = try acceptableShapeScore(coordinates.len);
    while (shapes.items.items.len != 0) {
        if (try matchShapes(allocator, shapes, data.constraints(), data.restraints(), &checked)) |matched| {
            found = true;
            if (matched.score > chosen_score) {
                if (chosen) |*shape| shape.deinit();
                chosen = try shapes.items.items[matched.shape].clone(allocator);
                chosen_start = matched.start;
                chosen_score = matched.score;
                if (matched.score > acceptable_score) break;
            }
        } else found = false;
        if (checked > max_macrocycles) break;
        var equivalents = try shapes.equivalents();
        errdefer equivalents.deinit();
        try equivalents.removeDuplicates();
        shapes.deinit();
        shapes = equivalents;
        equivalents = ShapeCollection.init(allocator);
    }
    if (found) if (chosen) |shape| {
        try writeCoordinates(shape, chosen_start, coordinates, coordinates_set);
        return if (shape.pentagon_vertices.items.len != 0)
            .matched_needs_minimization
        else
            .matched;
    };
    return .no_shape;
}

/// Topology-facing shape entry point used by the fragment layout layer. Atom
/// storage remains canonical; the temporary coordinate vectors follow the
/// ring path and are mapped back only after a successful match.
pub fn generateRingShape(
    allocator: std.mem.Allocator,
    ring: core.ids.RingId,
    atoms: []model.Atom,
    bonds: []const model.Bond,
    graph: topology.Graph,
    membership: topology.RingMembership,
    coordinates_set: []bool,
    force_open_macrocycles: bool,
) core.errors.Error!GenerateResult {
    if (coordinates_set.len != atoms.len) return error.InvalidMapping;
    var data = try collectPathData(allocator, ring, atoms, bonds, graph, membership);
    defer data.deinit();
    const coordinates = allocator.alloc(core.math.Vec2, data.ordered_atoms.len) catch return error.OutOfMemory;
    defer allocator.free(coordinates);
    const path_set = allocator.alloc(bool, data.ordered_atoms.len) catch return error.OutOfMemory;
    defer allocator.free(path_set);
    for (data.ordered_atoms, coordinates, path_set) |atom, *coordinate, *is_set| {
        coordinate.* = atoms[atom.index()].coordinates;
        is_set.* = coordinates_set[atom.index()];
    }
    const result = try generateShape(allocator, data, coordinates, path_set, force_open_macrocycles);
    if (result.placed()) for (data.ordered_atoms, coordinates) |atom, coordinate| {
        atoms[atom.index()].coordinates = coordinate;
        coordinates_set[atom.index()] = true;
    };
    return result;
}

/// Select the least disruptive bond for the recursive acyclic fallback. The
/// score and strict first-wins tie behavior mirror pinned upstream.
pub fn findBondToOpen(
    ring: core.ids.RingId,
    bonds: []const model.Bond,
    graph: topology.Graph,
    membership: topology.RingMembership,
) ?core.ids.BondId {
    if (!ring.isValid() or ring.index() >= membership.rings.len) return null;
    const is_macrocycle = membership.atoms(ring).len >= topology.rings.macrocycle_size;
    var best: ?core.ids.BondId = null;
    var best_score: f32 = 0;
    for (membership.ringBonds(ring)) |bond_id| {
        const bond = bonds[bond_id.index()];
        if (is_macrocycle) {
            if (bond.effective_order != .single) continue;
            var next_to_stereo = false;
            for (graph.incidentBonds(bond.start)) |other| if (bonds[other.index()].stereo != .unspecified) {
                next_to_stereo = true;
                break;
            };
            if (!next_to_stereo) for (graph.incidentBonds(bond.end)) |other| if (bonds[other.index()].stereo != .unspecified) {
                next_to_stereo = true;
                break;
            };
            if (next_to_stereo) continue;
        }
        if (!std.math.isFinite(bond.crossing_penalty_multiplier) or bond.crossing_penalty_multiplier <= 0) continue;
        const score = @as(f32, @floatFromInt(membership.bondRings(bond_id).len * 10 + graph.degree(bond.start) + graph.degree(bond.end))) /
            bond.crossing_penalty_multiplier;
        if (best == null or score < best_score) {
            best = bond_id;
            best_score = score;
        }
    }
    return best;
}

/// Emit the invert-bond degrees of freedom created for eligible macrocycle
/// substituents after neighbor placement. Affected atoms are the bound-side
/// submolecule in deterministic graph traversal order.
pub fn collectDofs(
    allocator: std.mem.Allocator,
    bonds: []const model.Bond,
    graph: topology.Graph,
    membership: topology.RingMembership,
    fragmentation: fragments.Fragmentation,
) core.errors.Error!core.dof.Collection {
    var dofs: std.ArrayList(core.dof.Dof) = .empty;
    defer dofs.deinit(allocator);
    var affected_atoms: std.ArrayList(core.ids.AtomId) = .empty;
    defer affected_atoms.deinit(allocator);
    for (fragmentation.atom_fragment, 0..) |fragment, atom_index| {
        if (!fragment.isValid()) continue;
        const atom = core.ids.AtomId.fromIndex(@intCast(atom_index));
        const atom_rings = membership.atomRings(atom);
        if (atom_rings.len != 1 or membership.atoms(atom_rings[0]).len < topology.rings.macrocycle_size or graph.degree(atom) != 3) continue;
        var blocked_by_stereo = false;
        for (graph.incidentBonds(atom)) |bond_id| {
            const bond = bonds[bond_id.index()];
            if (bond.stereo == .unspecified or graph.degree(bond.start) == 1 or graph.degree(bond.end) == 1) continue;
            blocked_by_stereo = true;
            break;
        }
        if (blocked_by_stereo) continue;
        for (graph.neighbors(atom)) |neighbor| {
            if (shareRing(atom, neighbor, membership)) continue;
            // Upstream's constructor adds exactly the bound atom. It does not
            // propagate through the substituent or across fragment frames.
            if (affected_atoms.items.len > std.math.maxInt(u32) or dofs.items.len >= std.math.maxInt(u32)) return error.TooManyItems;
            const affected_start: u32 = @intCast(affected_atoms.items.len);
            affected_atoms.append(allocator, neighbor) catch return error.OutOfMemory;
            dofs.append(allocator, .{
                .id = core.ids.DofId.fromIndex(@intCast(dofs.items.len)),
                .fragment = fragment,
                .affected_atoms = .{ .start = affected_start, .len = 1 },
                .state = .{ .count = 2, .tier = core.dof.Tier.invert_bond },
                .payload = .{ .invert_bond = .{ .pivot = atom, .bound = neighbor } },
            }) catch return error.OutOfMemory;
        }
    }
    const owned_dofs = dofs.toOwnedSlice(allocator) catch return error.OutOfMemory;
    errdefer allocator.free(owned_dofs);
    const owned_affected = affected_atoms.toOwnedSlice(allocator) catch return error.OutOfMemory;
    return .{ .allocator = allocator, .items = owned_dofs, .affected_atoms = owned_affected };
}

/// Collect the always-present primary fragment DOFs together with specialized
/// macrocycle DOFs in upstream fragment-local order.
pub fn collectAllDofs(
    allocator: std.mem.Allocator,
    atoms: []const model.Atom,
    bonds: []const model.Bond,
    graph: topology.Graph,
    membership: topology.RingMembership,
    fragmentation: fragments.Fragmentation,
) core.errors.Error!core.dof.Collection {
    var items: std.ArrayList(core.dof.Dof) = .empty;
    defer items.deinit(allocator);
    var affected: std.ArrayList(core.ids.AtomId) = .empty;
    defer affected.deinit(allocator);
    const maximum_specialized = std.math.mul(usize, bonds.len, 4) catch return error.TooManyItems;
    const primary_count = std.math.mul(usize, fragmentation.fragments.len, 3) catch return error.TooManyItems;
    const maximum_items = std.math.add(usize, primary_count, maximum_specialized) catch return error.TooManyItems;
    const maximum_affected = std.math.mul(usize, maximum_specialized, fragmentation.atom_fragment.len) catch return error.TooManyItems;
    items.ensureTotalCapacity(allocator, maximum_items) catch return error.OutOfMemory;
    affected.ensureTotalCapacity(allocator, maximum_affected) catch return error.OutOfMemory;

    for (fragmentation.fragments) |fragment| {
        const primary = [_]struct { payload: core.dof.Payload, count: u32, tier: u32 }{
            .{
                .payload = .{ .flip_fragment = .{} },
                .count = if (fragment.parent.isValid()) 2 else 1,
                .tier = core.dof.Tier.flip_fragment,
            },
            .{
                .payload = .{ .change_parent_bond_length = .{} },
                .count = 7,
                .tier = core.dof.Tier.change_parent_bond_length,
            },
            .{
                .payload = .{ .rotate_fragment = .{} },
                .count = if (fragment.parent.isValid()) 5 else 1,
                .tier = core.dof.Tier.rotate_fragment,
            },
        };
        for (primary) |entry| {
            if (items.items.len >= std.math.maxInt(u32) or affected.items.len > std.math.maxInt(u32)) return error.TooManyItems;
            items.append(allocator, .{
                .id = core.ids.DofId.fromIndex(@intCast(items.items.len)),
                .fragment = fragment.id,
                .affected_atoms = .{ .start = @intCast(affected.items.len) },
                .state = .{ .count = entry.count, .tier = entry.tier },
                .payload = entry.payload,
            }) catch return error.OutOfMemory;
        }
        try appendSpecializedDofs(allocator, &items, &affected, fragment, atoms, bonds, graph, membership, fragmentation);
    }
    const owned_items = items.toOwnedSlice(allocator) catch return error.OutOfMemory;
    errdefer allocator.free(owned_items);
    const owned_affected = affected.toOwnedSlice(allocator) catch return error.OutOfMemory;
    return .{ .allocator = allocator, .items = owned_items, .affected_atoms = owned_affected };
}

const PendingDof = struct {
    payload: core.dof.Payload,
    state: core.dof.State,
    affected: std.ArrayList(core.ids.AtomId) = .empty,
};

fn appendSpecializedDofs(
    allocator: std.mem.Allocator,
    output: *std.ArrayList(core.dof.Dof),
    output_affected: *std.ArrayList(core.ids.AtomId),
    fragment: fragments.Fragment,
    atoms: []const model.Atom,
    bonds: []const model.Bond,
    graph: topology.Graph,
    membership: topology.RingMembership,
    fragmentation: fragments.Fragmentation,
) core.errors.Error!void {
    var pending: std.ArrayList(PendingDof) = .empty;
    defer {
        for (pending.items) |*dof| dof.affected.deinit(allocator);
        pending.deinit(allocator);
    }
    const maximum_pending = std.math.mul(usize, bonds.len, 4) catch return error.TooManyItems;
    pending.ensureTotalCapacity(allocator, maximum_pending) catch return error.OutOfMemory;
    const membership_count = std.math.mul(usize, maximum_pending, fragmentation.atom_fragment.len) catch return error.TooManyItems;
    const atom_dofs = allocator.alloc(bool, membership_count) catch return error.OutOfMemory;
    defer allocator.free(atom_dofs);
    @memset(atom_dofs, false);
    const visited = allocator.alloc(bool, fragmentation.atom_fragment.len) catch return error.OutOfMemory;
    defer allocator.free(visited);
    @memset(visited, false);
    const queue = allocator.alloc(core.ids.AtomId, fragment.atom_count) catch return error.OutOfMemory;
    defer allocator.free(queue);
    const ordered_neighbors = allocator.alloc(core.ids.AtomId, fragmentation.atom_fragment.len) catch return error.OutOfMemory;
    defer allocator.free(ordered_neighbors);
    var head: usize = 0;
    var tail: usize = 0;
    const members = fragmentation.members(fragment.id);
    for (members) |atom| {
        if (membership.atomRings(atom).len == 0) continue;
        queue[tail] = atom;
        tail += 1;
        visited[atom.index()] = true;
    }
    if (tail == 0) {
        var start = members[0];
        if (fragment.parent.isValid()) {
            // The attachment end, resolved by Fragmentation rather than read
            // off the bond's stored direction (cgz-jg4).
            if (fragment.attachment_atom.isValid()) start = fragment.attachment_atom;
        }
        queue[0] = start;
        tail = 1;
        visited[start.index()] = true;
    }

    while (head < tail) : (head += 1) {
        const atom = queue[head];
        const neighbors = graph.neighbors(atom);
        if (membership.atomRings(atom).len == 0) {
            try neighbour_order.orderNeighbours(allocator, atoms, bonds, graph, membership, atom, ordered_neighbors[0..neighbors.len]);
            var start: usize = 0;
            for (ordered_neighbors[0..neighbors.len], 0..) |neighbor, index| {
                if (visited[neighbor.index()]) {
                    start = index;
                    break;
                }
            }
            std.mem.rotate(core.ids.AtomId, ordered_neighbors[0..neighbors.len], start);
        } else {
            @memcpy(ordered_neighbors[0..neighbors.len], neighbors);
        }
        for (ordered_neighbors[0..neighbors.len]) |neighbor| {
            if (visited[neighbor.index()]) continue;
            if (fragmentation.atom_fragment[neighbor.index()] == fragment.id) {
                visited[neighbor.index()] = true;
                queue[tail] = neighbor;
                tail += 1;
            }
            // Upstream adds every newly placed neighbour to the current
            // atom's existing DOFs before deciding whether that neighbour is
            // part of this fragment's traversal. A child-fragment attachment
            // therefore belongs to the affected set even though it is not
            // queued here.
            for (0..pending.items.len) |pending_index| {
                if (!atom_dofs[pending_index * fragmentation.atom_fragment.len + atom.index()]) continue;
                pending.items[pending_index].affected.append(allocator, neighbor) catch return error.OutOfMemory;
                atom_dofs[pending_index * fragmentation.atom_fragment.len + neighbor.index()] = true;
            }
        }

        if (eligibleInvertPivot(atom, bonds, graph, membership)) {
            for (graph.neighbors(atom)) |neighbor| {
                if (shareRing(atom, neighbor, membership)) continue;
                try appendPendingDof(allocator, &pending, atom_dofs, fragmentation.atom_fragment.len, fragment.atom_count, .{ .invert_bond = .{ .pivot = atom, .bound = neighbor } }, .{
                    .count = 2,
                    .tier = core.dof.Tier.invert_bond,
                }, neighbor);
            }
        }
        for (graph.neighbors(atom)) |neighbor| {
            if (fragmentation.atom_fragment[neighbor.index()] != fragment.id or shareRing(atom, neighbor, membership)) continue;
            try appendPendingDof(allocator, &pending, atom_dofs, fragmentation.atom_fragment.len, fragment.atom_count, .{ .scale_atoms = .{ .pivot = atom } }, .{
                .count = 2,
                .tier = core.dof.Tier.scale_atoms,
            }, neighbor);
        }
    }
    if (tail != members.len) return error.InvalidMapping;
    for (pending.items) |dof| {
        if (output.items.len >= std.math.maxInt(u32) or output_affected.items.len > std.math.maxInt(u32) or dof.affected.items.len > std.math.maxInt(u32)) return error.TooManyItems;
        const affected_start: u32 = @intCast(output_affected.items.len);
        output_affected.appendSlice(allocator, dof.affected.items) catch return error.OutOfMemory;
        output.append(allocator, .{
            .id = core.ids.DofId.fromIndex(@intCast(output.items.len)),
            .fragment = fragment.id,
            .affected_atoms = .{ .start = affected_start, .len = @intCast(dof.affected.items.len) },
            .state = dof.state,
            .payload = dof.payload,
        }) catch return error.OutOfMemory;
    }
}

fn appendPendingDof(
    allocator: std.mem.Allocator,
    pending: *std.ArrayList(PendingDof),
    atom_dofs: []bool,
    atom_count: usize,
    maximum_affected: usize,
    payload: core.dof.Payload,
    state: core.dof.State,
    first_atom: core.ids.AtomId,
) core.errors.Error!void {
    const index = pending.items.len;
    pending.append(allocator, .{ .payload = payload, .state = state }) catch return error.OutOfMemory;
    pending.items[index].affected.ensureTotalCapacity(allocator, maximum_affected) catch return error.OutOfMemory;
    pending.items[index].affected.append(allocator, first_atom) catch return error.OutOfMemory;
    atom_dofs[index * atom_count + first_atom.index()] = true;
}

fn eligibleInvertPivot(atom: core.ids.AtomId, bonds: []const model.Bond, graph: topology.Graph, membership: topology.RingMembership) bool {
    const atom_rings = membership.atomRings(atom);
    if (atom_rings.len != 1 or membership.atoms(atom_rings[0]).len < topology.rings.macrocycle_size or graph.degree(atom) != 3) return false;
    for (graph.incidentBonds(atom)) |bond_id| {
        const bond = bonds[bond_id.index()];
        if (bond.stereo != .unspecified and graph.degree(bond.start) != 1 and graph.degree(bond.end) != 1) return false;
    }
    return true;
}

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

fn bondBetween(left: core.ids.AtomId, right: core.ids.AtomId, graph: topology.Graph) ?core.ids.BondId {
    for (graph.neighbors(left), graph.incidentBonds(left)) |neighbor, bond| {
        if (neighbor == right) return bond;
    }
    return null;
}

fn shareRing(left: core.ids.AtomId, right: core.ids.AtomId, membership: topology.RingMembership) bool {
    for (membership.atomRings(left)) |ring| {
        if (std.mem.indexOfScalar(core.ids.RingId, membership.atomRings(right), ring) != null) return true;
    }
    return false;
}

fn orderRing(
    allocator: std.mem.Allocator,
    ring: core.ids.RingId,
    graph: topology.Graph,
    membership: topology.RingMembership,
) core.errors.Error![]core.ids.AtomId {
    if (!ring.isValid() or ring.index() >= membership.rings.len) return error.InvalidMapping;
    const members = membership.atoms(ring);
    if (members.len < 3) return error.InvalidMapping;
    const ordered = allocator.alloc(core.ids.AtomId, members.len) catch return error.OutOfMemory;
    errdefer allocator.free(ordered);
    ordered[0] = members[0];
    var previous = core.ids.AtomId.invalid;
    var current = members[0];
    for (1..members.len) |index| {
        var next: ?core.ids.AtomId = null;
        for (graph.neighbors(current), graph.incidentBonds(current)) |neighbor, bond_id| {
            if (neighbor == previous or std.mem.indexOfScalar(core.ids.RingId, membership.bondRings(bond_id), ring) == null) continue;
            next = neighbor;
            break;
        }
        const selected = next orelse return error.InvalidMapping;
        if (selected == ordered[0] or std.mem.indexOfScalar(core.ids.AtomId, ordered[0..index], selected) != null) return error.InvalidMapping;
        ordered[index] = selected;
        previous = current;
        current = selected;
    }
    const closing_bond = bondBetween(current, ordered[0], graph) orelse return error.InvalidMapping;
    if (std.mem.indexOfScalar(core.ids.RingId, membership.bondRings(closing_bond), ring) == null) return error.InvalidMapping;
    return ordered;
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

test "round polyomino supports the nine-atom macrocycle boundary" {
    var polyomino = Polyomino.init(std.testing.allocator);
    defer polyomino.deinit();
    try polyomino.buildWithVertices(9);
    const contour = try polyomino.path(std.testing.allocator);
    defer std.testing.allocator.free(contour);
    try std.testing.expectEqual(@as(usize, 9), contour.len);
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
    try core.oom.checkAllocationFailures(std.testing.allocator, buildAndTraverse, .{});
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
    try core.oom.checkAllocationFailures(std.testing.allocator, enumerateAndDeduplicate, .{});
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
    try core.oom.checkAllocationFailures(std.testing.allocator, matchAndWrite, .{});
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

fn collectAndGenerateFixture(allocator: std.mem.Allocator, test_dofs: bool) !void {
    const ring_size = 13;
    var atoms: [ring_size + 2]model.Atom = undefined;
    for (&atoms, 0..) |*atom, index| atom.* = .{
        .id = core.ids.AtomId.fromIndex(@intCast(index)),
        .input_index = @intCast(index),
        .atomic_number = if (index == 4) .nitrogen else .carbon,
    };
    var bonds: [ring_size + 2]model.Bond = undefined;
    for (0..ring_size) |index| bonds[index] = .{
        .id = core.ids.BondId.fromIndex(@intCast(index)),
        .input_index = @intCast(index),
        .start = core.ids.AtomId.fromIndex(@intCast(index)),
        .end = core.ids.AtomId.fromIndex(@intCast((index + 1) % ring_size)),
        .input_order = .single,
        .effective_order = .single,
    };
    bonds[0].input_order = .double;
    bonds[0].effective_order = .double;
    bonds[0].stereo = .z;
    bonds[0].stereo_atom_a = core.ids.AtomId.fromIndex(ring_size - 1);
    bonds[0].stereo_atom_b = core.ids.AtomId.fromIndex(2);
    bonds[ring_size] = .{
        .id = core.ids.BondId.fromIndex(ring_size),
        .input_index = ring_size,
        .start = core.ids.AtomId.fromIndex(2),
        .end = core.ids.AtomId.fromIndex(ring_size),
        .input_order = .single,
        .effective_order = .single,
    };
    bonds[ring_size + 1] = .{
        .id = core.ids.BondId.fromIndex(ring_size + 1),
        .input_index = ring_size + 1,
        .start = core.ids.AtomId.fromIndex(ring_size),
        .end = core.ids.AtomId.fromIndex(ring_size + 1),
        .input_order = .single,
        .effective_order = .single,
    };
    var graph = try topology.Graph.init(allocator, &atoms, &bonds);
    defer graph.deinit();
    var membership = try topology.RingMembership.init(allocator, graph, &bonds);
    defer membership.deinit();
    try std.testing.expectEqual(@as(usize, 1), membership.rings.len);
    const opened = findBondToOpen(core.ids.RingId.fromIndex(0), &bonds, graph, membership) orelse return error.InvalidMapping;
    try std.testing.expect(bonds[opened.index()].start != core.ids.AtomId.fromIndex(0));
    try std.testing.expect(bonds[opened.index()].end != core.ids.AtomId.fromIndex(0));
    try std.testing.expect(bonds[opened.index()].start != core.ids.AtomId.fromIndex(1));
    try std.testing.expect(bonds[opened.index()].end != core.ids.AtomId.fromIndex(1));
    if (test_dofs) {
        var fragmentation = try fragments.Fragmentation.init(allocator, &atoms, &bonds, graph, membership);
        defer fragmentation.deinit();
        var dofs = try collectDofs(allocator, &bonds, graph, membership, fragmentation);
        defer dofs.deinit();
        try std.testing.expectEqual(@as(usize, 1), dofs.items.len);
        try std.testing.expectEqual(core.dof.DofKind.invert_bond, dofs.items[0].kind());
        try std.testing.expectEqual(core.ids.AtomId.fromIndex(2), dofs.items[0].payload.invert_bond.pivot);
        try std.testing.expectEqual(core.ids.AtomId.fromIndex(ring_size), dofs.items[0].payload.invert_bond.bound);
        try std.testing.expectEqualSlices(core.ids.AtomId, &.{
            core.ids.AtomId.fromIndex(ring_size),
        }, dofs.affected_atoms);
        try core.oom.checkAllocationFailures(
            std.testing.allocator,
            collectDofsAndDiscard,
            .{ &bonds, graph, membership, fragmentation },
        );
        var all_dofs = try collectAllDofs(allocator, &atoms, &bonds, graph, membership, fragmentation);
        defer all_dofs.deinit();
        try std.testing.expect(all_dofs.items.len > fragmentation.fragments.len * 3 + 1);
        var found_scale = false;
        var found_invert = false;
        var dof_index: usize = 0;
        for (fragmentation.fragments) |fragment| {
            try std.testing.expectEqual(core.dof.DofKind.flip_fragment, all_dofs.items[dof_index].kind());
            try std.testing.expectEqual(core.dof.DofKind.change_parent_bond_length, all_dofs.items[dof_index + 1].kind());
            try std.testing.expectEqual(core.dof.DofKind.rotate_fragment, all_dofs.items[dof_index + 2].kind());
            try std.testing.expectEqual(fragment.id, all_dofs.items[dof_index].fragment);
            dof_index += 3;
            while (dof_index < all_dofs.items.len and all_dofs.items[dof_index].fragment == fragment.id) : (dof_index += 1) {
                found_scale = found_scale or all_dofs.items[dof_index].kind() == .scale_atoms;
                found_invert = found_invert or all_dofs.items[dof_index].kind() == .invert_bond;
            }
        }
        try std.testing.expect(found_scale);
        try std.testing.expect(found_invert);
        try core.oom.checkAllocationFailures(
            std.testing.allocator,
            collectAllDofsAndDiscard,
            .{ &atoms, &bonds, graph, membership, fragmentation },
        );
    }
    var data = try collectPathData(allocator, core.ids.RingId.fromIndex(0), &atoms, &bonds, graph, membership);
    defer data.deinit();
    try std.testing.expectEqual(@as(usize, ring_size), data.ordered_atoms.len);
    try std.testing.expectEqualSlices(usize, &.{4}, data.hetero_atoms);
    try std.testing.expectEqual(@as(usize, 1), data.substituted_atoms.len);
    try std.testing.expectEqual(@as(usize, 2), data.substituted_atoms[0].subtree_size);
    try std.testing.expectEqual(@as(usize, 1), data.double_bonds.len);
    try std.testing.expect(!data.double_bonds[0].trans);

    const coordinates = try allocator.alloc(core.math.Vec2, ring_size);
    defer allocator.free(coordinates);
    @memset(coordinates, .{});
    const coordinates_set = try allocator.alloc(bool, ring_size);
    defer allocator.free(coordinates_set);
    @memset(coordinates_set, false);
    // Asserts a shape was applied. The pentagon-vertex distinction is a
    // separate fact reported by the same result, and this fixture does not
    // pin which of the two it lands on.
    try std.testing.expect((try generateShape(allocator, data, coordinates, coordinates_set, false)).placed());
    for (coordinates, 0..) |coordinate, index| {
        const following = coordinates[(index + 1) % coordinates.len];
        const dx = coordinate.x - following.x;
        const dy = coordinate.y - following.y;
        // Hex contour edges are √3 * BONDLENGTH. An odd ring omits one
        // contour vertex to model a pentagon, producing two BONDLENGTH edges.
        const length = @sqrt(dx * dx + dy * dy);
        const hex_edge = core.math.bond_length * @sqrt(@as(f32, 3));
        try std.testing.expect(@abs(length - core.math.bond_length) < 0.001 or @abs(length - hex_edge) < 0.001);
    }
    const before = try allocator.dupe(core.math.Vec2, coordinates);
    defer allocator.free(before);
    try std.testing.expectEqual(GenerateResult.no_shape, try generateShape(allocator, data, coordinates, coordinates_set, true));
    try std.testing.expectEqualSlices(core.math.Vec2, before, coordinates);
}

fn collectDofsAndDiscard(
    allocator: std.mem.Allocator,
    bonds: []const model.Bond,
    graph: topology.Graph,
    membership: topology.RingMembership,
    fragmentation: fragments.Fragmentation,
) !void {
    var dofs = try collectDofs(allocator, bonds, graph, membership, fragmentation);
    defer dofs.deinit();
}

fn collectAllDofsAndDiscard(
    allocator: std.mem.Allocator,
    atoms: []const model.Atom,
    bonds: []const model.Bond,
    graph: topology.Graph,
    membership: topology.RingMembership,
    fragmentation: fragments.Fragmentation,
) !void {
    var dofs = try collectAllDofs(allocator, atoms, bonds, graph, membership, fragmentation);
    defer dofs.deinit();
}

test "five-atom path DOFs match the stable upstream probe" {
    const input_atoms = [_]topology.prepare.TestAtom{ .{}, .{}, .{}, .{}, .{} };
    const input_bonds = [_]topology.prepare.TestBond{
        .{ .start = 0, .end = 1 },
        .{ .start = 1, .end = 2 },
        .{ .start = 2, .end = 3 },
        .{ .start = 3, .end = 4 },
    };
    var prepared = try topology.prepareInput(std.testing.allocator, topology.prepare.TestInput{
        .atoms = &input_atoms,
        .bonds = &input_bonds,
    });
    defer prepared.deinit();
    var fragmentation = try fragments.Fragmentation.init(
        std.testing.allocator,
        prepared.working.atoms,
        prepared.working.bonds,
        prepared.graph,
        prepared.rings,
    );
    defer fragmentation.deinit();
    var dofs = try collectAllDofs(
        std.testing.allocator,
        prepared.working.atoms,
        prepared.working.bonds,
        prepared.graph,
        prepared.rings,
        fragmentation,
    );
    defer dofs.deinit();

    const expected_kinds = [_]core.dof.DofKind{
        .flip_fragment,   .change_parent_bond_length, .rotate_fragment, .scale_atoms,   .scale_atoms,
        .flip_fragment,   .change_parent_bond_length, .rotate_fragment, .flip_fragment, .change_parent_bond_length,
        .rotate_fragment, .scale_atoms,               .scale_atoms,
    };
    const expected_fragments = [_]u32{ 0, 0, 0, 0, 0, 1, 1, 1, 2, 2, 2, 2, 2 };
    const expected_counts = [_]u32{ 1, 7, 1, 2, 2, 2, 7, 5, 2, 7, 5, 2, 2 };
    const expected_tiers = [_]u32{ 0, 2, 3, 4, 4, 0, 2, 3, 0, 2, 3, 4, 4 };
    const expected_affected = [_]u32{ 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1 };
    const expected_affected_starts = [_]u32{ 0, 0, 0, 0, 1, 2, 2, 2, 2, 2, 2, 2, 3 };
    const expected_affected_atoms = [_]u32{ 0, 1, 4, 3 };
    const expected_pivots = [_]u32{
        std.math.maxInt(u32), std.math.maxInt(u32), std.math.maxInt(u32), 1,                    0,
        std.math.maxInt(u32), std.math.maxInt(u32), std.math.maxInt(u32), std.math.maxInt(u32), std.math.maxInt(u32),
        std.math.maxInt(u32), 3,                    4,
    };
    try std.testing.expectEqual(expected_kinds.len, dofs.items.len);
    for (dofs.items, 0..) |dof, index| {
        try std.testing.expectEqual(core.ids.DofId.fromIndex(@intCast(index)), dof.id);
        try std.testing.expectEqual(expected_kinds[index], dof.kind());
        try std.testing.expectEqual(expected_fragments[index], dof.fragment.index());
        try std.testing.expectEqual(@as(u32, 0), dof.state.current);
        try std.testing.expectEqual(@as(u32, 0), dof.state.optimal);
        try std.testing.expectEqual(expected_counts[index], dof.state.count);
        try std.testing.expectEqual(expected_tiers[index], dof.state.tier);
        try std.testing.expectEqual(expected_affected_starts[index], dof.affected_atoms.start);
        try std.testing.expectEqual(expected_affected[index], dof.affected_atoms.len);
        const pivot = switch (dof.payload) {
            .scale_atoms => |scale| prepared.working.atoms[scale.pivot.index()].input_index,
            else => std.math.maxInt(u32),
        };
        try std.testing.expectEqual(expected_pivots[index], pivot);
    }
    for (dofs.affected_atoms, expected_affected_atoms) |atom, expected| {
        try std.testing.expectEqual(expected, prepared.working.atoms[atom.index()].input_index);
    }
}

test "topology path extraction and bounded shape orchestration clean every allocation failure" {
    // ArrayList growth can remap in place through the debug allocator on
    // macOS depending on the address reused by the preceding injected-failure
    // run. Page allocation makes those growth decisions stable while
    // FailingAllocator still verifies every allocation index and byte balance.
    try core.oom.checkAllocationFailures(std.heap.page_allocator, collectAndGenerateFixture, .{false});
}

test "macrocycle substituents emit invert-bond DOFs for exactly the bound atom" {
    try collectAndGenerateFixture(std.testing.allocator, true);
}
