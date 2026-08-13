const std = @import("std");
const core = @import("core");
const model = @import("model");
const data = @import("templates/data.zig");

pub const template_count = data.templates.len;

pub const Edge = struct {
    start: u8,
    end: u8,
    order: core.chemistry.BondOrder = .single,
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
        if (try compare(allocator, atom_count, bonds, templateEdges(template_index))) |mapping| {
            return .{
                .allocator = allocator,
                .template_index = @intCast(template_index),
                .mapping = mapping,
            };
        }
    }
    return null;
}

/// Apply normalized immutable coordinates at the conserved 50-unit bond
/// length. Candidate atoms are supplied in the same order used by findGraph.
pub fn apply(match: Match, candidate_atoms: []const core.ids.AtomId, atoms: []model.Atom) core.errors.Error!void {
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

fn compare(
    allocator: std.mem.Allocator,
    atom_count: usize,
    molecule_bonds: []const Edge,
    template_bonds: []const data.Bond,
) core.errors.Error!?[]u8 {
    const molecule_scores = allocator.alloc(u32, atom_count) catch return error.OutOfMemory;
    defer allocator.free(molecule_scores);
    const template_scores = allocator.alloc(u32, atom_count) catch return error.OutOfMemory;
    defer allocator.free(template_scores);
    if (try morganScores(allocator, atom_count, molecule_bonds, molecule_scores) !=
        try morganScores(allocator, atom_count, template_bonds, template_scores)) return null;

    const matrix = allocator.alloc(bool, atom_count * atom_count) catch return error.OutOfMemory;
    defer allocator.free(matrix);
    for (molecule_scores, 0..) |molecule_score, molecule_index| {
        for (template_scores, 0..) |template_score, template_index| {
            matrix[molecule_index * atom_count + template_index] = molecule_score == template_score;
        }
    }
    const mapping = allocator.alloc(u8, atom_count) catch return error.OutOfMemory;
    errdefer allocator.free(mapping);
    const used = allocator.alloc(bool, atom_count) catch return error.OutOfMemory;
    defer allocator.free(used);
    @memset(used, false);
    if (!searchIdentity(0, mapping, used, matrix, atom_count, molecule_bonds, template_bonds)) {
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
    template_bonds: []const data.Bond,
) bool {
    if (depth == atom_count) return true;
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
        if (searchIdentity(depth + 1, mapping, used, matrix, atom_count, molecule_bonds, template_bonds)) return true;
        used[template_atom] = false;
    }
    return false;
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
    try apply(match, &ids, &atoms);
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
