const std = @import("std");
const core = @import("core");

pub const Atom = struct {
    atomic_number: core.chemistry.AtomicNumber,
    hidden: bool = false,
};

pub const Bond = struct {
    start: core.ids.InputIndex,
    end: core.ids.InputIndex,
    effective_order: core.chemistry.BondOrder,
    skip: bool = false,
};

pub const Result = struct {
    allocator: std.mem.Allocator,
    /// Visible atoms in upstream canonical order, followed by hidden atoms in
    /// caller order. This is a complete permutation suitable for owned model
    /// storage; `visible_count` delimits the upstream-observable prefix.
    internal_to_input: []core.ids.InputIndex,
    scores: []u32,
    structural_bonds: []core.ids.InputIndex,
    visible_count: u32,

    pub fn deinit(self: *Result) void {
        self.allocator.free(self.structural_bonds);
        self.allocator.free(self.scores);
        self.allocator.free(self.internal_to_input);
        self.* = undefined;
    }
};

/// Port of sketcherMinimizer::morganScores and canonicalOrdering. Stable tie
/// handling and the upstream cumulative `newScores` buffer are intentional.
pub fn order(
    allocator: std.mem.Allocator,
    atoms: []const Atom,
    bonds: []const Bond,
) core.errors.Error!Result {
    if (atoms.len > std.math.maxInt(u32) or bonds.len > std.math.maxInt(u32)) {
        return error.TooManyItems;
    }

    const input_to_visible = allocator.alloc(u32, atoms.len) catch return error.OutOfMemory;
    defer allocator.free(input_to_visible);
    @memset(input_to_visible, std.math.maxInt(u32));
    var visible_count: usize = 0;
    for (atoms, 0..) |atom, input_index| {
        if (atom.hidden) continue;
        input_to_visible[input_index] = @intCast(visible_count);
        visible_count += 1;
    }

    const old_scores = allocator.alloc(u32, visible_count) catch return error.OutOfMemory;
    defer allocator.free(old_scores);
    @memset(old_scores, 1);
    if (visible_count < 2) @memset(old_scores, 0);
    if (visible_count >= 2) try morganScores(allocator, bonds, input_to_visible, old_scores);

    const canonical_scores = allocator.alloc(u32, visible_count) catch return error.OutOfMemory;
    defer allocator.free(canonical_scores);
    for (canonical_scores, 0..) |*score, index| {
        const input_index = visibleInputAt(input_to_visible, @intCast(index));
        score.* = old_scores[index] *% 100 +% @backingInt(atoms[input_index].atomic_number);
    }

    const permutation = allocator.alloc(core.ids.InputIndex, atoms.len) catch return error.OutOfMemory;
    errdefer allocator.free(permutation);
    const output_scores = allocator.alloc(u32, atoms.len) catch return error.OutOfMemory;
    errdefer allocator.free(output_scores);
    @memset(output_scores, 0);
    var structural_count: usize = 0;
    for (bonds) |bond| structural_count += @intFromBool(isStructural(bond, input_to_visible));
    const structural_bonds = allocator.alloc(core.ids.InputIndex, structural_count) catch return error.OutOfMemory;
    errdefer allocator.free(structural_bonds);

    const atom_visited = allocator.alloc(bool, visible_count) catch return error.OutOfMemory;
    defer allocator.free(atom_visited);
    @memset(atom_visited, false);
    const bond_visited = allocator.alloc(bool, bonds.len) catch return error.OutOfMemory;
    defer allocator.free(bond_visited);
    @memset(bond_visited, false);
    const queue = allocator.alloc(core.ids.InputIndex, visible_count) catch return error.OutOfMemory;
    defer allocator.free(queue);

    var output_count: usize = 0;
    var bond_output_count: usize = 0;
    while (output_count < visible_count) {
        var root: ?u32 = null;
        for (0..visible_count) |candidate| {
            if (atom_visited[candidate]) continue;
            if (root == null or canonical_scores[candidate] > canonical_scores[root.?]) {
                root = @intCast(candidate);
            }
        }
        var head: usize = 0;
        var tail: usize = 1;
        queue[0] = visibleInputAt(input_to_visible, root.?);
        atom_visited[root.?] = true;
        while (head < tail) : (head += 1) {
            const current = queue[head];
            const current_visible = input_to_visible[current];
            permutation[output_count] = current;
            output_scores[output_count] = old_scores[current_visible];
            output_count += 1;

            while (bestUnvisitedBond(bonds, bond_visited, input_to_visible, canonical_scores, current)) |bond_index| {
                bond_visited[bond_index] = true;
                structural_bonds[bond_output_count] = @intCast(bond_index);
                bond_output_count += 1;
                const bond = bonds[bond_index];
                const neighbor = if (bond.start == current) bond.end else bond.start;
                const neighbor_visible = input_to_visible[neighbor];
                if (!atom_visited[neighbor_visible]) {
                    atom_visited[neighbor_visible] = true;
                    queue[tail] = neighbor;
                    tail += 1;
                }
            }
        }
    }

    for (atoms, 0..) |atom, input_index| {
        if (!atom.hidden) continue;
        permutation[output_count] = @intCast(input_index);
        output_count += 1;
    }
    std.debug.assert(output_count == atoms.len);
    std.debug.assert(bond_output_count == structural_bonds.len);
    return .{
        .allocator = allocator,
        .internal_to_input = permutation,
        .scores = output_scores,
        .structural_bonds = structural_bonds,
        .visible_count = @intCast(visible_count),
    };
}

fn morganScores(
    allocator: std.mem.Allocator,
    bonds: []const Bond,
    input_to_visible: []const u32,
    old_scores: []u32,
) core.errors.Error!void {
    const new_scores = allocator.alloc(u32, old_scores.len) catch return error.OutOfMemory;
    defer allocator.free(new_scores);
    @memset(new_scores, 0);
    const ordered_scores = allocator.alloc(u32, old_scores.len) catch return error.OutOfMemory;
    defer allocator.free(ordered_scores);
    var old_ties = old_scores.len;
    while (true) {
        for (bonds) |bond| {
            if (!isStructural(bond, input_to_visible)) continue;
            const start = input_to_visible[bond.start];
            const end = input_to_visible[bond.end];
            new_scores[start] +%= old_scores[end];
            new_scores[end] +%= old_scores[start];
        }
        @memcpy(ordered_scores, new_scores);
        std.mem.sort(u32, ordered_scores, {}, std.sort.asc(u32));
        var new_ties: usize = 0;
        for (ordered_scores[1..], ordered_scores[0 .. ordered_scores.len - 1]) |score, previous| {
            new_ties += @intFromBool(score == previous);
        }
        if (new_ties >= old_ties) break;
        old_ties = new_ties;
        @memcpy(old_scores, new_scores);
    }
}

fn bestUnvisitedBond(
    bonds: []const Bond,
    visited: []const bool,
    input_to_visible: []const u32,
    scores: []const u32,
    current: core.ids.InputIndex,
) ?usize {
    var best: ?usize = null;
    for (bonds, 0..) |bond, index| {
        if (visited[index] or !isStructural(bond, input_to_visible)) continue;
        if (bond.start != current and bond.end != current) continue;
        const neighbor = if (bond.start == current) bond.end else bond.start;
        if (best) |best_index| {
            const best_bond = bonds[best_index];
            const best_neighbor = if (best_bond.start == current) best_bond.end else best_bond.start;
            if (scores[input_to_visible[neighbor]] <= scores[input_to_visible[best_neighbor]]) continue;
        }
        best = index;
    }
    return best;
}

fn isStructural(bond: Bond, input_to_visible: []const u32) bool {
    return !bond.skip and bond.effective_order != .zero and
        bond.start < input_to_visible.len and bond.end < input_to_visible.len and
        input_to_visible[bond.start] != std.math.maxInt(u32) and
        input_to_visible[bond.end] != std.math.maxInt(u32);
}

fn visibleInputAt(input_to_visible: []const u32, visible_index: u32) core.ids.InputIndex {
    for (input_to_visible, 0..) |candidate, input_index| {
        if (candidate == visible_index) return @intCast(input_index);
    }
    unreachable;
}

test "Morgan ordering is deterministic and hidden atoms follow the visible prefix" {
    const atoms = [_]Atom{
        .{ .atomic_number = .carbon },
        .{ .atomic_number = .oxygen },
        .{ .atomic_number = .carbon, .hidden = true },
        .{ .atomic_number = .nitrogen },
    };
    const bonds = [_]Bond{
        .{ .start = 0, .end = 1, .effective_order = .single },
        .{ .start = 0, .end = 3, .effective_order = .single },
        .{ .start = 1, .end = 2, .effective_order = .single },
    };
    var result = try order(std.testing.allocator, &atoms, &bonds);
    defer result.deinit();
    try std.testing.expectEqual(@as(u32, 3), result.visible_count);
    try std.testing.expectEqualSlices(u32, &.{ 0, 1, 3, 2 }, result.internal_to_input);
    try std.testing.expectEqualSlices(u32, &.{ 2, 1, 1, 0 }, result.scores);
}

test "canonical cycle and disconnected orders match pinned oracle probes" {
    const cycle_atoms = [_]Atom{
        .{ .atomic_number = .carbon },
        .{ .atomic_number = .nitrogen },
        .{ .atomic_number = .oxygen },
        .{ .atomic_number = .carbon },
    };
    const cycle_bonds = [_]Bond{
        .{ .start = 0, .end = 1, .effective_order = .single },
        .{ .start = 1, .end = 2, .effective_order = .single },
        .{ .start = 2, .end = 3, .effective_order = .single },
        .{ .start = 3, .end = 0, .effective_order = .single },
    };
    var cycle = try order(std.testing.allocator, &cycle_atoms, &cycle_bonds);
    defer cycle.deinit();
    try std.testing.expectEqualSlices(u32, &.{ 2, 1, 3, 0 }, cycle.internal_to_input);
    try std.testing.expectEqualSlices(u32, &.{ 2, 2, 2, 2 }, cycle.scores);

    const disconnected_atoms = [_]Atom{
        .{ .atomic_number = .carbon },
        .{ .atomic_number = .oxygen },
        .{ .atomic_number = .nitrogen },
        .{ .atomic_number = .carbon },
        .{ .atomic_number = .sulfur },
    };
    const disconnected_bonds = [_]Bond{
        .{ .start = 0, .end = 1, .effective_order = .single },
        .{ .start = 2, .end = 3, .effective_order = .single },
    };
    var disconnected = try order(std.testing.allocator, &disconnected_atoms, &disconnected_bonds);
    defer disconnected.deinit();
    try std.testing.expectEqualSlices(u32, &.{ 1, 0, 2, 3, 4 }, disconnected.internal_to_input);
}

test "Morgan ordering reports and cleans every allocation failure" {
    const atoms = [_]Atom{ .{ .atomic_number = .carbon }, .{ .atomic_number = .oxygen } };
    const bonds = [_]Bond{.{ .start = 0, .end = 1, .effective_order = .single }};
    const Runner = struct {
        fn run(allocator: std.mem.Allocator, atom_items: []const Atom, bond_items: []const Bond) !void {
            var result = try order(allocator, atom_items, bond_items);
            result.deinit();
        }
    };
    try core.oom.checkAllocationFailures(std.testing.allocator, Runner.run, .{ &atoms, &bonds });
}
