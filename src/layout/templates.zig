const std = @import("std");
const core = @import("core");
const model = @import("model");
const topology = @import("topology");
/// The pinned template fixture. Public because the regeneration check
/// (tests/template_generate.zig) validates these exact tables, and giving
/// that file its own module over the same source would place one file in two
/// modules — which Zig rejects as soon as anything reaches the layout layer.
pub const data = @import("templates/data.zig");

pub const template_count = data.templates.len;

pub const Options = struct {
    load_templates: bool = true,
    /// Runtime MAE templates are unavailable in the dependency-free build and
    /// therefore rejected explicitly, matching the stable oracle adapter.
    template_directory: ?[]const u8 = null,
};

pub const Edge = struct {
    start: u8,
    end: u8,
    order: core.chemistry.BondOrder = .single,
    /// Absolute Z/E requirement. Callers leave this unspecified for ordinary
    /// and small-ring double bonds, matching upstream's ignore-Z/E rule.
    stereo: core.chemistry.BondStereo = .unspecified,
    start_neighbor: ?u8 = null,
    end_neighbor: ?u8 = null,
};

/// One immutable-database match. `mapping[candidate_atom]` is the template
/// atom index, exactly matching the upstream probe's direction.
pub const Match = struct {
    allocator: std.mem.Allocator,
    template_index: u8,
    mapping: []u8,

    pub fn deinit(self: *Match) void {
        self.allocator.free(self.mapping);
        self.* = undefined;
    }
};

/// A match plus the native ring-system atom order used to obtain it. Both
/// slices are context-owned; no scratch is stored in the immutable database.
pub const RingMatch = struct {
    allocator: std.mem.Allocator,
    atoms: []core.ids.AtomId,
    match: Match,

    pub fn deinit(self: *RingMatch) void {
        self.match.deinit();
        self.allocator.free(self.atoms);
        self.* = undefined;
    }

    pub fn apply(self: RingMatch, working_atoms: []model.Atom) core.errors.Error!void {
        try applyMatch(self.match, self.atoms, working_atoms);
    }

    pub fn applyAligned(self: RingMatch, working_atoms: []model.Atom) core.errors.Error!void {
        try self.apply(working_atoms);
        alignConstrained(working_atoms, self.atoms);
    }
};

pub fn templateEdges(index: usize) []const data.Bond {
    const template = data.templates[index];
    return data.bonds[template.bond_start..][0..template.bond_len];
}

pub fn templateCoordinates(index: usize, atom_index: usize) core.math.Vec2 {
    const template = data.templates[index];
    const atom = data.atoms[template.atom_start + atom_index];
    return .{ .x = @bitCast(atom.x_bits), .y = @bitCast(atom.y_bits) };
}

/// Find the first embedded template using upstream's cumulative Morgan scores
/// and depth-first identity traversal. The embedded database is compile-time
/// immutable; all score, matrix, and traversal state is allocator-owned by the
/// calling generation context.
pub fn findGraph(
    allocator: std.mem.Allocator,
    atom_count: usize,
    bonds: []const Edge,
) core.errors.Error!?Match {
    if (atom_count > std.math.maxInt(u8)) return error.TooManyItems;
    for (bonds) |bond| {
        if (bond.start >= atom_count or bond.end >= atom_count or bond.start == bond.end) {
            return error.InvalidAtomIndex;
        }
    }
    for (data.templates, 0..) |template, template_index| {
        if (template.atom_len != atom_count) continue;
        if (try compare(allocator, atom_count, bonds, template_index)) |mapping| {
            return .{
                .allocator = allocator,
                .template_index = @intCast(template_index),
                .mapping = mapping,
            };
        }
    }
    return null;
}

/// Match one upstream-style ring set. Ring and atom traversal order is
/// retained exactly while shared fused atoms and bonds are emitted once.
/// `load_templates = false` is a real no-op rather than a lazy-init switch.
pub fn findRingSet(
    allocator: std.mem.Allocator,
    options: Options,
    atom_count: usize,
    bonds: []const model.Bond,
    membership: anytype,
    ring_ids: []const core.ids.RingId,
) core.errors.Error!?RingMatch {
    if (options.template_directory != null) return error.Unsupported;
    if (!options.load_templates or ring_ids.len < 2) return null;
    const local_index = allocator.alloc(u8, atom_count) catch return error.OutOfMemory;
    defer allocator.free(local_index);
    @memset(local_index, std.math.maxInt(u8));
    var ring_atoms: std.ArrayList(core.ids.AtomId) = .empty;
    defer ring_atoms.deinit(allocator);
    for (ring_ids) |ring_id| {
        for (membership.atoms(ring_id)) |atom_id| {
            if (atom_id.index() >= atom_count) return error.InvalidAtomIndex;
            if (local_index[atom_id.index()] != std.math.maxInt(u8)) continue;
            if (ring_atoms.items.len >= std.math.maxInt(u8)) return error.TooManyItems;
            local_index[atom_id.index()] = @intCast(ring_atoms.items.len);
            ring_atoms.append(allocator, atom_id) catch return error.OutOfMemory;
        }
    }

    const seen_bonds = allocator.alloc(bool, bonds.len) catch return error.OutOfMemory;
    defer allocator.free(seen_bonds);
    @memset(seen_bonds, false);
    var edges: std.ArrayList(Edge) = .empty;
    defer edges.deinit(allocator);
    for (ring_ids) |ring_id| {
        for (membership.ringBonds(ring_id)) |bond_id| {
            if (bond_id.index() >= bonds.len) return error.InvalidMapping;
            if (seen_bonds[bond_id.index()]) continue;
            seen_bonds[bond_id.index()] = true;
            const bond = bonds[bond_id.index()];
            const start = local_index[bond.start.index()];
            const end = local_index[bond.end.index()];
            if (start == std.math.maxInt(u8) or end == std.math.maxInt(u8)) return error.InvalidMapping;
            edges.append(allocator, .{ .start = start, .end = end, .order = bond.effective_order }) catch return error.OutOfMemory;
        }
    }
    var match = (try findGraph(allocator, ring_atoms.items.len, edges.items)) orelse return null;
    errdefer match.deinit();
    const owned_atoms = ring_atoms.toOwnedSlice(allocator) catch return error.OutOfMemory;
    return .{ .allocator = allocator, .atoms = owned_atoms, .match = match };
}

/// Apply normalized immutable coordinates at the conserved 50-unit bond
/// length. Candidate atoms are supplied in the same order used by findGraph.
pub fn applyMatch(match: Match, candidate_atoms: []const core.ids.AtomId, atoms: []model.Atom) core.errors.Error!void {
    if (candidate_atoms.len != match.mapping.len) return error.InvalidMapping;
    const template = data.templates[match.template_index];
    for (candidate_atoms, match.mapping) |atom_id, template_atom| {
        if (atom_id.index() >= atoms.len or template_atom >= template.atom_len) return error.InvalidMapping;
        const coordinate = templateCoordinates(match.template_index, template_atom);
        atoms[atom_id.index()].coordinates = .{
            .x = coordinate.x * core.math.bond_length,
            .y = coordinate.y * core.math.bond_length,
        };
    }
}

/// Rigid least-squares 2D alignment against caller constraints, followed by
/// exact restoration of fixed coordinates. This is allocation-free and keeps
/// constraint targets immutable.
pub fn alignConstrained(atoms: []model.Atom, candidate_atoms: []const core.ids.AtomId) void {
    var source_center: core.math.Vec2 = .{};
    var target_center: core.math.Vec2 = .{};
    var count: usize = 0;
    for (candidate_atoms) |atom_id| {
        const atom = atoms[atom_id.index()];
        if (!atom.constrained and !atom.fixed) continue;
        const target = atom.template_coordinates orelse continue;
        source_center.x += atom.coordinates.x;
        source_center.y += atom.coordinates.y;
        target_center.x += target.x;
        target_center.y += target.y;
        count += 1;
    }
    if (count != 0) {
        const inverse_count = 1 / @as(f32, @floatFromInt(count));
        source_center.x *= inverse_count;
        source_center.y *= inverse_count;
        target_center.x *= inverse_count;
        target_center.y *= inverse_count;
        var dot: f32 = 0;
        var cross: f32 = 0;
        for (candidate_atoms) |atom_id| {
            const atom = atoms[atom_id.index()];
            if (!atom.constrained and !atom.fixed) continue;
            const target = atom.template_coordinates orelse continue;
            const sx = atom.coordinates.x - source_center.x;
            const sy = atom.coordinates.y - source_center.y;
            const tx = target.x - target_center.x;
            const ty = target.y - target_center.y;
            dot += sx * tx + sy * ty;
            cross += sx * ty - sy * tx;
        }
        const angle = if (count > 1) std.math.atan2(cross, dot) else 0;
        const cosine = @cos(angle);
        const sine = @sin(angle);
        for (candidate_atoms) |atom_id| {
            const atom = &atoms[atom_id.index()];
            const x = atom.coordinates.x - source_center.x;
            const y = atom.coordinates.y - source_center.y;
            atom.coordinates = .{
                .x = x * cosine - y * sine + target_center.x,
                .y = x * sine + y * cosine + target_center.y,
            };
        }
    }
    for (candidate_atoms) |atom_id| {
        const atom = &atoms[atom_id.index()];
        if (atom.fixed) {
            if (atom.template_coordinates) |target| atom.coordinates = target;
        }
    }
}

pub fn rmsd(reference: []const core.math.Vec2, points: []const core.math.Vec2) core.errors.Error!f32 {
    if (reference.len != points.len) return error.InvalidMapping;
    var total: f32 = 0;
    for (reference, points) |expected, actual| {
        const dx = expected.x - actual.x;
        const dy = expected.y - actual.y;
        total += dx * dx + dy * dy;
    }
    if (reference.len != 0) total /= @floatFromInt(reference.len);
    return @sqrt(total);
}

/// Upstream's rounded 2x2 SVD alignment matrix. Keeping this exact seam is
/// important because template constraint orientation observes its rounding.
pub fn alignmentMatrix(reference: []const core.math.Vec2, points: []const core.math.Vec2) core.errors.Error![4]f32 {
    if (reference.len != points.len) return error.InvalidMapping;
    var a = [4]f32{ 0, 0, 0, 0 };
    for (reference, points) |ref, point| {
        a[0] += ref.x * point.x;
        a[1] += ref.y * point.x;
        a[2] += ref.x * point.y;
        a[3] += ref.y * point.y;
    }
    const decomposition = svd2x2(a);
    const u = decomposition.u;
    const v = decomposition.v;
    return .{
        v[0] * u[0] + v[1] * u[1],
        v[0] * u[2] + v[1] * u[3],
        v[2] * u[0] + v[3] * u[1],
        v[2] * u[2] + v[3] * u[3],
    };
}

const Svd = struct { u: [4]f32, singular: [4]f32, v: [4]f32 };

fn svd2x2(a: [4]f32) Svd {
    const transpose = [4]f32{ a[0], a[2], a[1], a[3] };
    const su = [4]f32{
        a[0] * transpose[0] + a[1] * transpose[2],
        a[0] * transpose[1] + a[1] * transpose[3],
        a[2] * transpose[0] + a[3] * transpose[2],
        a[2] * transpose[1] + a[3] * transpose[3],
    };
    const phi: f32 = 0.5 * std.math.atan2(su[1] + su[2], su[0] - su[3]);
    const cphi = roundTwo(@cos(phi));
    const sphi = roundTwo(@sin(phi));
    const u = [4]f32{ -cphi, -sphi, -sphi, cphi };

    const sw = [4]f32{
        transpose[0] * a[0] + transpose[1] * a[2],
        transpose[0] * a[1] + transpose[1] * a[3],
        transpose[2] * a[0] + transpose[3] * a[2],
        transpose[2] * a[1] + transpose[3] * a[3],
    };
    const theta: f32 = 0.5 * std.math.atan2(sw[1] + sw[2], sw[0] - sw[3]);
    const ctheta = @cos(theta);
    const stheta = @sin(theta);
    const w = [4]f32{ ctheta, -stheta, stheta, ctheta };
    const sum = su[0] + su[3];
    const difference = @sqrt((su[0] - su[3]) * (su[0] - su[3]) + 4 * su[1] * su[2]);
    const singular = [4]f32{ @sqrt((sum + difference) * 0.5), 0, 0, @sqrt((sum - difference) * 0.5) };

    const u_transpose = [4]f32{ u[0], u[2], u[1], u[3] };
    const ua = [4]f32{
        u_transpose[0] * a[0] + u_transpose[1] * a[2],
        u_transpose[0] * a[1] + u_transpose[1] * a[3],
        u_transpose[2] * a[0] + u_transpose[3] * a[2],
        u_transpose[2] * a[1] + u_transpose[3] * a[3],
    };
    const s = [4]f32{
        roundTwo(ua[0] * w[0] + ua[1] * w[2]),
        roundTwo(ua[0] * w[1] + ua[1] * w[3]),
        roundTwo(ua[2] * w[0] + ua[3] * w[2]),
        roundTwo(ua[2] * w[1] + ua[3] * w[3]),
    };
    const sign1: f32 = if (s[0] < 0) -1 else 1;
    const sign2: f32 = if (s[3] < 0) -1 else 1;
    return .{
        .u = u,
        .singular = singular,
        .v = .{
            roundTwo(w[0] * sign1),
            roundTwo(w[1] * sign2),
            roundTwo(w[2] * sign1),
            roundTwo(w[3] * sign2),
        },
    };
}

fn roundTwo(value: f32) f32 {
    return @round(value * 100) / 100;
}

fn compare(
    allocator: std.mem.Allocator,
    atom_count: usize,
    molecule_bonds: []const Edge,
    template_index: usize,
) core.errors.Error!?[]u8 {
    const template_bonds = templateEdges(template_index);
    const molecule_scores = allocator.alloc(u32, atom_count) catch return error.OutOfMemory;
    defer allocator.free(molecule_scores);
    const template_scores = allocator.alloc(u32, atom_count) catch return error.OutOfMemory;
    defer allocator.free(template_scores);
    if (try morganScores(allocator, atom_count, molecule_bonds, molecule_scores) !=
        try morganScores(allocator, atom_count, template_bonds, template_scores)) return null;

    const matrix = allocator.alloc(bool, atom_count * atom_count) catch return error.OutOfMemory;
    defer allocator.free(matrix);
    for (molecule_scores, 0..) |molecule_score, molecule_index| {
        for (template_scores, 0..) |template_score, template_atom| {
            matrix[molecule_index * atom_count + template_atom] = molecule_score == template_score;
        }
    }
    const mapping = allocator.alloc(u8, atom_count) catch return error.OutOfMemory;
    errdefer allocator.free(mapping);
    const used = allocator.alloc(bool, atom_count) catch return error.OutOfMemory;
    defer allocator.free(used);
    @memset(used, false);
    if (!searchIdentity(0, mapping, used, matrix, atom_count, molecule_bonds, template_index)) {
        allocator.free(mapping);
        return null;
    }
    return mapping;
}

/// This deliberately retains upstream's cumulative new-score buffer and tie
/// stopping rule rather than substituting a conventional Morgan refinement.
fn morganScores(
    allocator: std.mem.Allocator,
    atom_count: usize,
    bonds: anytype,
    old_scores: []u32,
) core.errors.Error!u32 {
    if (atom_count < 2) return 0;
    @memset(old_scores, 1);
    const new_scores = allocator.alloc(u32, atom_count) catch return error.OutOfMemory;
    defer allocator.free(new_scores);
    @memset(new_scores, 0);
    const ordered = allocator.alloc(u32, atom_count) catch return error.OutOfMemory;
    defer allocator.free(ordered);
    var old_ties = atom_count;
    var iterations: u32 = 0;
    while (true) {
        iterations += 1;
        for (bonds) |bond| {
            const start = edgeStart(bond);
            const end = edgeEnd(bond);
            new_scores[start] +%= old_scores[end];
            new_scores[end] +%= old_scores[start];
        }
        @memcpy(ordered, new_scores);
        std.mem.sort(u32, ordered, {}, std.sort.asc(u32));
        var new_ties: usize = 0;
        for (ordered[1..], ordered[0 .. ordered.len - 1]) |score, previous| {
            new_ties += @intFromBool(score == previous);
        }
        if (new_ties >= old_ties) return iterations;
        old_ties = new_ties;
        @memcpy(old_scores, new_scores);
    }
}

fn edgeStart(bond: anytype) u8 {
    return if (@hasField(@TypeOf(bond), "start")) bond.start else bond.from;
}

fn edgeEnd(bond: anytype) u8 {
    return if (@hasField(@TypeOf(bond), "end")) bond.end else bond.to;
}

fn searchIdentity(
    depth: usize,
    mapping: []u8,
    used: []bool,
    matrix: []const bool,
    atom_count: usize,
    molecule_bonds: []const Edge,
    template_index: usize,
) bool {
    const template_bonds = templateEdges(template_index);
    if (depth == atom_count) return stereoMappingMatches(mapping, molecule_bonds, template_index);
    for (0..atom_count) |template_atom| {
        if (!matrix[depth * atom_count + template_atom] or used[template_atom]) continue;
        var compatible = true;
        for (molecule_bonds) |bond| {
            const high = @max(bond.start, bond.end);
            if (high != depth) continue;
            const low = @min(bond.start, bond.end);
            if (!hasTemplateBond(template_bonds, @intCast(template_atom), mapping[low])) {
                compatible = false;
                break;
            }
        }
        if (!compatible) continue;
        mapping[depth] = @intCast(template_atom);
        used[template_atom] = true;
        if (searchIdentity(depth + 1, mapping, used, matrix, atom_count, molecule_bonds, template_index)) return true;
        used[template_atom] = false;
    }
    return false;
}

fn stereoMappingMatches(mapping: []const u8, bonds: []const Edge, template_index: usize) bool {
    for (bonds) |bond| {
        const expected_z = switch (bond.stereo) {
            .cis, .z => true,
            .trans, .e => false,
            .unspecified => continue,
        };
        const start_neighbor = bond.start_neighbor orelse continue;
        const end_neighbor = bond.end_neighbor orelse continue;
        if (start_neighbor >= mapping.len or end_neighbor >= mapping.len) return false;
        const p1 = templateCoordinates(template_index, mapping[start_neighbor]);
        const p2 = templateCoordinates(template_index, mapping[bond.start]);
        const p3 = templateCoordinates(template_index, mapping[bond.end]);
        const p4 = templateCoordinates(template_index, mapping[end_neighbor]);
        const line_x = p3.x - p2.x;
        const line_y = p3.y - p2.y;
        const side1 = line_x * (p1.y - p2.y) - line_y * (p1.x - p2.x);
        const side4 = line_x * (p4.y - p2.y) - line_y * (p4.x - p2.x);
        if ((side1 * side4 > 0) != expected_z) return false;
    }
    return true;
}

fn hasTemplateBond(bonds: []const data.Bond, left: u8, right: u8) bool {
    for (bonds) |bond| {
        if ((bond.from == left and bond.to == right) or (bond.from == right and bond.to == left)) return true;
    }
    return false;
}

fn fixtureEdges(comptime template_index: usize) [data.templates[template_index].bond_len]Edge {
    const bonds = templateEdges(template_index);
    var result: [bonds.len]Edge = undefined;
    for (bonds, 0..) |bond, index| result[index] = .{
        .start = bond.from,
        .end = bond.to,
        .order = core.chemistry.BondOrder.fromInt(bond.order).?,
    };
    return result;
}

fn findAndDiscard(allocator: std.mem.Allocator, edges: []const Edge) !void {
    var match = (try findGraph(allocator, 5, edges)) orelse return error.Unsupported;
    match.deinit();
}

test "five-atom fused graph selects pinned template 80 with exact mapping" {
    const edges = comptime fixtureEdges(80);
    var match = (try findGraph(std.testing.allocator, 5, &edges)).?;
    defer match.deinit();
    try std.testing.expectEqual(@as(u8, 80), match.template_index);
    try std.testing.expectEqualSlices(u8, &.{ 0, 1, 2, 3, 4 }, match.mapping);

    var atoms: [5]model.Atom = undefined;
    var ids: [5]core.ids.AtomId = undefined;
    for (&atoms, &ids, 0..) |*atom, *id, index| {
        id.* = core.ids.AtomId.fromIndex(@intCast(index));
        atom.* = .{ .id = id.*, .input_index = @intCast(index), .atomic_number = .carbon };
    }
    try applyMatch(match, &ids, &atoms);
    for (atoms, match.mapping) |atom, template_atom| {
        const expected = templateCoordinates(80, template_atom);
        try std.testing.expectEqual(expected.x * core.math.bond_length, atom.coordinates.x);
        try std.testing.expectEqual(expected.y * core.math.bond_length, atom.coordinates.y);
    }
}

test "template matching rejects nonmatching and malformed graphs" {
    const path = [_]Edge{ .{ .start = 0, .end = 1 }, .{ .start = 1, .end = 2 } };
    try std.testing.expect((try findGraph(std.testing.allocator, 3, &path)) == null);
    try std.testing.expectError(error.InvalidAtomIndex, findGraph(
        std.testing.allocator,
        2,
        &.{.{ .start = 0, .end = 2 }},
    ));
}

test "template matching cleans every injected allocation failure" {
    const edges = comptime fixtureEdges(80);
    try std.testing.checkAllAllocationFailures(std.testing.allocator, findAndDiscard, .{&edges});
}

test "all 82 immutable templates are structurally matchable" {
    var edge_buffer: [64]Edge = undefined;
    for (data.templates, 0..) |template, template_index| {
        const source_edges = templateEdges(template_index);
        try std.testing.expect(source_edges.len <= edge_buffer.len);
        for (source_edges, edge_buffer[0..source_edges.len]) |source, *destination| {
            destination.* = .{
                .start = source.from,
                .end = source.to,
                .order = core.chemistry.BondOrder.fromInt(source.order).?,
            };
        }
        var match = (try findGraph(std.testing.allocator, template.atom_len, edge_buffer[0..source_edges.len])).?;
        defer match.deinit();
        // Earlier entries intentionally win when the database contains an
        // isomorphic graph, exactly as in upstream's linear search.
        try std.testing.expect(match.template_index <= template_index);
        var seen: [256]bool = @splat(false);
        for (match.mapping) |template_atom| {
            try std.testing.expect(template_atom < data.templates[match.template_index].atom_len);
            try std.testing.expect(!seen[template_atom]);
            seen[template_atom] = true;
        }
        for (edge_buffer[0..source_edges.len]) |edge| {
            try std.testing.expect(hasTemplateBond(
                templateEdges(match.template_index),
                match.mapping[edge.start],
                match.mapping[edge.end],
            ));
        }
    }
}

test "template identity enforces absolute double-bond geometry" {
    const mapping = [_]u8{ 0, 1, 2, 3, 4 };
    var edge = Edge{
        .start = 0,
        .end = 1,
        .order = .double,
        .stereo = .z,
        .start_neighbor = 3,
        .end_neighbor = 2,
    };
    const identity_is_z = stereoMappingMatches(&mapping, &.{edge}, 80);
    edge.stereo = .e;
    try std.testing.expect(identity_is_z != stereoMappingMatches(&mapping, &.{edge}, 80));
}

fn fixtureAtom(index: u32) model.Atom {
    return .{ .id = core.ids.AtomId.fromIndex(index), .input_index = index, .atomic_number = .carbon };
}

fn fixtureBond(index: u32, start: u32, end: u32) model.Bond {
    return .{
        .id = core.ids.BondId.fromIndex(index),
        .input_index = index,
        .start = core.ids.AtomId.fromIndex(start),
        .end = core.ids.AtomId.fromIndex(end),
        .input_order = .single,
        .effective_order = .single,
    };
}

fn prepareTemplate80(
    allocator: std.mem.Allocator,
    atoms: *[5]model.Atom,
    bonds: *[6]model.Bond,
) !struct { topology.Graph, topology.RingMembership } {
    for (atoms, 0..) |*atom, index| atom.* = fixtureAtom(@intCast(index));
    const edges = comptime fixtureEdges(80);
    for (bonds, edges, 0..) |*bond, edge, index| bond.* = fixtureBond(@intCast(index), edge.start, edge.end);
    var graph = try topology.Graph.init(allocator, atoms, bonds);
    errdefer graph.deinit();
    const rings = try topology.RingMembership.init(allocator, graph, bonds);
    return .{ graph, rings };
}

test "native fused ring set selects and applies embedded template 80" {
    var atoms: [5]model.Atom = undefined;
    var bonds: [6]model.Bond = undefined;
    var graph, var rings = try prepareTemplate80(std.testing.allocator, &atoms, &bonds);
    defer graph.deinit();
    defer rings.deinit();
    var ring_ids: [2]core.ids.RingId = undefined;
    try std.testing.expectEqual(ring_ids.len, rings.rings.len);
    for (&ring_ids, 0..) |*ring_id, index| ring_id.* = core.ids.RingId.fromIndex(@intCast(index));
    var match = (try findRingSet(std.testing.allocator, .{}, atoms.len, &bonds, rings, &ring_ids)).?;
    defer match.deinit();
    try std.testing.expectEqual(@as(u8, 80), match.match.template_index);
    try match.apply(&atoms);
    for (match.atoms, match.match.mapping) |atom_id, template_atom| {
        const expected = templateCoordinates(80, template_atom);
        try std.testing.expectEqual(expected.x * core.math.bond_length, atoms[atom_id.index()].coordinates.x);
        try std.testing.expectEqual(expected.y * core.math.bond_length, atoms[atom_id.index()].coordinates.y);
    }
    try std.testing.expect((try findRingSet(std.testing.allocator, .{ .load_templates = false }, atoms.len, &bonds, rings, &ring_ids)) == null);
    try std.testing.expectError(error.Unsupported, findRingSet(
        std.testing.allocator,
        .{ .template_directory = "user-templates" },
        atoms.len,
        &bonds,
        rings,
        &ring_ids,
    ));
}

test "template application aligns constraints and restores fixed coordinates" {
    const edges = comptime fixtureEdges(80);
    var match = (try findGraph(std.testing.allocator, 5, &edges)).?;
    defer match.deinit();
    var atoms: [5]model.Atom = undefined;
    var ids: [5]core.ids.AtomId = undefined;
    for (&atoms, &ids, 0..) |*atom, *id, index| {
        id.* = core.ids.AtomId.fromIndex(@intCast(index));
        atom.* = fixtureAtom(@intCast(index));
    }
    const offset = core.math.Vec2{ .x = 17, .y = -23 };
    for (0..2) |index| {
        const source = templateCoordinates(80, match.mapping[index]);
        atoms[index].constrained = true;
        atoms[index].template_coordinates = .{
            .x = -source.y * core.math.bond_length + offset.x,
            .y = source.x * core.math.bond_length + offset.y,
        };
    }
    atoms[0].fixed = true;
    try applyMatch(match, &ids, &atoms);
    alignConstrained(&atoms, &ids);
    for (atoms, match.mapping) |atom, template_atom| {
        const source = templateCoordinates(80, template_atom);
        try std.testing.expectApproxEqAbs(-source.y * core.math.bond_length + offset.x, atom.coordinates.x, 0.001);
        try std.testing.expectApproxEqAbs(source.x * core.math.bond_length + offset.y, atom.coordinates.y, 0.001);
    }
    try std.testing.expectEqual(atoms[0].template_coordinates.?, atoms[0].coordinates);
}

test "RMSD and rounded SVD alignment retain upstream f32 behavior" {
    const reference = [_]core.math.Vec2{ .{ .x = -1 }, .{ .x = 1 } };
    const rotated = [_]core.math.Vec2{ .{ .y = -1 }, .{ .y = 1 } };
    try std.testing.expectEqual(@as(f32, @sqrt(2.0)), try rmsd(&reference, &rotated));
    const matrix = try alignmentMatrix(&reference, &rotated);
    for (matrix) |value| try std.testing.expect(std.math.isFinite(value));
    const x = rotated[1].x * matrix[0] + rotated[1].y * matrix[1];
    const y = rotated[1].x * matrix[2] + rotated[1].y * matrix[3];
    try std.testing.expectApproxEqAbs(reference[1].x, x, 0.01);
    try std.testing.expectApproxEqAbs(reference[1].y, y, 0.01);
    try std.testing.expectError(error.InvalidMapping, rmsd(&reference, rotated[0..1]));
}

fn findRingSetAndDiscard(
    allocator: std.mem.Allocator,
    bonds: []const model.Bond,
    rings: topology.RingMembership,
    ring_ids: []const core.ids.RingId,
) !void {
    var match = (try findRingSet(allocator, .{}, 5, bonds, rings, ring_ids)) orelse return error.Unsupported;
    match.deinit();
}

test "ring-set matching cleans every injected allocation failure" {
    var atoms: [5]model.Atom = undefined;
    var bonds: [6]model.Bond = undefined;
    var graph, var rings = try prepareTemplate80(std.testing.allocator, &atoms, &bonds);
    defer graph.deinit();
    defer rings.deinit();
    const ring_ids = [_]core.ids.RingId{ core.ids.RingId.fromIndex(0), core.ids.RingId.fromIndex(1) };
    try std.testing.checkAllAllocationFailures(std.testing.allocator, findRingSetAndDiscard, .{ &bonds, rings, &ring_ids });
}

const ThreadResult = struct {
    err: ?anyerror = null,
    coordinates: [5]core.math.Vec2 = undefined,
};

fn renderTemplate80(result: *ThreadResult) void {
    var buffer: [64 * 1024]u8 = undefined;
    var fixed = std.heap.FixedBufferAllocator.init(&buffer);
    const allocator = fixed.allocator();
    var atoms: [5]model.Atom = undefined;
    var bonds: [6]model.Bond = undefined;
    var graph, var rings = prepareTemplate80(allocator, &atoms, &bonds) catch |err| {
        result.err = err;
        return;
    };
    defer graph.deinit();
    defer rings.deinit();
    const ring_ids = [_]core.ids.RingId{ core.ids.RingId.fromIndex(0), core.ids.RingId.fromIndex(1) };
    var match = (findRingSet(allocator, .{}, atoms.len, &bonds, rings, &ring_ids) catch |err| {
        result.err = err;
        return;
    }) orelse {
        result.err = error.Unsupported;
        return;
    };
    defer match.deinit();
    match.apply(&atoms) catch |err| {
        result.err = err;
        return;
    };
    for (atoms) |atom| result.coordinates[atom.input_index] = atom.coordinates;
}

test "two first-use template contexts produce serial-identical coordinates" {
    var serial: ThreadResult = .{};
    var first: ThreadResult = .{};
    var second: ThreadResult = .{};
    renderTemplate80(&serial);
    const first_thread = try std.Thread.spawn(.{}, renderTemplate80, .{&first});
    const second_thread = try std.Thread.spawn(.{}, renderTemplate80, .{&second});
    first_thread.join();
    second_thread.join();
    try std.testing.expect(serial.err == null);
    try std.testing.expect(first.err == null);
    try std.testing.expect(second.err == null);
    try std.testing.expectEqual(serial.coordinates, first.coordinates);
    try std.testing.expectEqual(serial.coordinates, second.coordinates);
}
