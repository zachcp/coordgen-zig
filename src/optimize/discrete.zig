const std = @import("std");
const core = @import("core");
const geometry = @import("geometry");

pub const PenaltyContext = struct {
    constrained_flip: bool = false,
    chain_with_chain_parent: bool = false,
};

pub const clash_energy_threshold: f32 = 10;
pub const rejected_solution_score: f32 = @floatFromInt(99_999_999);
pub const sketcher_epsilon: f32 = 0.0001;
pub const maximum_scored_solutions: usize = 10_000;

pub const Evaluator = struct {
    context: ?*anyopaque = null,
    scoreFn: *const fn (?*anyopaque, []const u32) anyerror!f32,

    fn score(self: Evaluator, states: []const u32) !f32 {
        return self.scoreFn(self.context, states);
    }
};

pub const SolutionCache = struct {
    allocator: std.mem.Allocator,
    state_count: usize,
    max_solutions: usize,
    evaluator: Evaluator,
    states: std.ArrayList(u32) = .empty,
    scores: std.ArrayList(f32) = .empty,

    pub fn init(allocator: std.mem.Allocator, state_count: usize, precision: f32, evaluator: Evaluator) core.errors.Error!SolutionCache {
        if (!std.math.isFinite(precision) or precision <= 0) return error.InvalidOption;
        const scaled = @as(f64, maximum_scored_solutions) * @as(f64, precision);
        if (scaled > @as(f64, @floatFromInt(std.math.maxInt(usize)))) return error.TooManyItems;
        return .{
            .allocator = allocator,
            .state_count = state_count,
            .max_solutions = @intFromFloat(scaled),
            .evaluator = evaluator,
        };
    }

    pub fn deinit(self: *SolutionCache) void {
        self.scores.deinit(self.allocator);
        self.states.deinit(self.allocator);
        self.* = undefined;
    }

    pub fn count(self: SolutionCache) usize {
        return self.scores.items.len;
    }

    pub fn score(self: *SolutionCache, candidate: []const u32) !f32 {
        if (candidate.len != self.state_count) return error.InvalidMapping;
        for (0..self.count()) |index| {
            if (std.mem.eql(u32, self.solution(index), candidate)) return self.scores.items[index];
        }
        // Preserve upstream's strict `>` guard: one entry beyond the scaled
        // nominal maximum is accepted before later unseen states are rejected.
        if (self.count() > self.max_solutions) return rejected_solution_score;
        const result = try self.evaluator.score(candidate);
        try self.states.ensureUnusedCapacity(self.allocator, candidate.len);
        try self.scores.ensureUnusedCapacity(self.allocator, 1);
        self.states.appendSliceAssumeCapacity(candidate);
        self.scores.appendAssumeCapacity(result);
        return result;
    }

    pub fn best(self: SolutionCache) core.errors.Error!struct { states: []const u32, score: f32 } {
        if (self.count() == 0) return error.InvalidMapping;
        var best_index: usize = 0;
        for (1..self.count()) |index| {
            const score_value = self.scores.items[index];
            if (score_value < self.scores.items[best_index] or
                (score_value == self.scores.items[best_index] and lexicographicLess(self.solution(index), self.solution(best_index))))
            {
                best_index = index;
            }
        }
        return .{ .states = self.solution(best_index), .score = self.scores.items[best_index] };
    }

    fn solution(self: SolutionCache, index: usize) []const u32 {
        return self.states.items[index * self.state_count ..][0..self.state_count];
    }
};

pub const ExhaustiveResult = struct {
    score: f32,
    clean_pose: bool,
};

pub fn exhaustiveSearch(
    allocator: std.mem.Allocator,
    dofs: []core.dof.Dof,
    selected: []const usize,
    cache: *SolutionCache,
    initial_score: f32,
) !ExhaustiveResult {
    if (cache.state_count != dofs.len) return error.InvalidMapping;
    for (selected) |index| if (index >= dofs.len) return error.InvalidMapping;
    const states = allocator.alloc(u32, dofs.len) catch return error.OutOfMemory;
    defer allocator.free(states);
    var context = ExhaustiveContext{
        .dofs = dofs,
        .selected = selected,
        .cache = cache,
        .states = states,
        .best_score = initial_score,
    };
    try context.visit(0);
    for (selected) |index| dofs[index].state.restoreOptimal();
    return .{ .score = context.best_score, .clean_pose = context.best_score < clash_energy_threshold };
}

const ExhaustiveContext = struct {
    dofs: []core.dof.Dof,
    selected: []const usize,
    cache: *SolutionCache,
    states: []u32,
    best_score: f32,
    abort: bool = false,

    fn visit(self: *ExhaustiveContext, depth: usize) !void {
        if (self.abort) return;
        if (depth == self.selected.len) {
            for (self.dofs, self.states) |dof, *state| state.* = dof.state.current;
            const result = try self.cache.score(self.states);
            if (result < clash_energy_threshold) {
                for (self.selected) |index| self.dofs[index].state.storeOptimal();
                self.best_score = result;
                self.abort = true;
            } else if (result < self.best_score - sketcher_epsilon) {
                self.best_score = result;
                for (self.selected) |index| self.dofs[index].state.storeOptimal();
            }
            return;
        }
        const index = self.selected[depth];
        const count = self.dofs[index].state.count;
        for (0..count) |_| {
            try self.visit(depth + 1);
            self.dofs[index].state.advance();
        }
    }
};

pub const Tuples = struct {
    allocator: std.mem.Allocator,
    width: usize,
    items: []usize,

    pub fn deinit(self: *Tuples) void {
        self.allocator.free(self.items);
        self.* = undefined;
    }

    pub fn count(self: Tuples) usize {
        return if (self.width == 0) 0 else self.items.len / self.width;
    }

    pub fn tuple(self: Tuples, index: usize) []const usize {
        return self.items[index * self.width ..][0..self.width];
    }
};

pub fn buildTuples(allocator: std.mem.Allocator, dof_count: usize, order: usize) core.errors.Error!Tuples {
    if (order == 0 or order > dof_count) return .{ .allocator = allocator, .width = order, .items = try allocator.alloc(usize, 0) };
    var output: std.ArrayList(usize) = .empty;
    defer output.deinit(allocator);
    const working = allocator.alloc(usize, order) catch return error.OutOfMemory;
    defer allocator.free(working);
    try appendTuples(allocator, &output, working, 0, 0, dof_count);
    return .{ .allocator = allocator, .width = order, .items = output.toOwnedSlice(allocator) catch return error.OutOfMemory };
}

fn appendTuples(allocator: std.mem.Allocator, output: *std.ArrayList(usize), working: []usize, depth: usize, start: usize, count: usize) core.errors.Error!void {
    if (depth == working.len) {
        output.appendSlice(allocator, working) catch return error.OutOfMemory;
        return;
    }
    const remaining = working.len - depth;
    for (start..count - remaining + 1) |index| {
        working[depth] = index;
        try appendTuples(allocator, output, working, depth + 1, index + 1, count);
    }
}

fn lexicographicLess(left: []const u32, right: []const u32) bool {
    for (left, right) |left_value, right_value| {
        if (left_value != right_value) return left_value < right_value;
    }
    return false;
}

pub fn affectedAtoms(collection: core.dof.Collection, dof: core.dof.Dof) core.errors.Error![]const core.ids.AtomId {
    const start = dof.affected_atoms.start;
    const end = std.math.add(u32, start, dof.affected_atoms.len) catch return error.InvalidMapping;
    if (end > collection.affected_atoms.len) return error.InvalidMapping;
    return collection.affected_atoms[start..end];
}

/// Apply one DOF to fragment-local coordinates. Callers rebuild global
/// fragment placement after applying all current states, as pinned upstream
/// does in CoordgenDOFSolutions::scoreCurrentSolution.
pub fn apply(
    dof: core.dof.Dof,
    affected: []const core.ids.AtomId,
    coordinates: []core.math.Vec2,
) !void {
    try dof.state.validate();
    for (affected) |atom| if (atom.index() >= coordinates.len) return error.InvalidAtomIndex;
    if (dof.state.current == 0) return;
    switch (dof.payload) {
        .flip_fragment => {
            for (affected) |atom| coordinates[atom.index()].y = -coordinates[atom.index()].y;
        },
        .change_parent_bond_length => {
            const factor = stateFactor(1.6, dof.state.current);
            const move_by = core.math.bond_length * (factor - 1);
            for (affected) |atom| coordinates[atom.index()].x += move_by;
        },
        .rotate_fragment => {
            var angle = @as(f32, std.math.pi) / 180 * 15 * @as(f32, @floatFromInt((dof.state.current + 1) / 2));
            if (dof.state.current % 2 == 0) angle = -angle;
            const origin = core.math.Vec2{ .x = -core.math.bond_length };
            for (affected) |atom| {
                const offset = subtract(coordinates[atom.index()], origin);
                coordinates[atom.index()] = add(geometry.rotate(offset, @sin(angle), @cos(angle)), origin);
            }
        },
        .scale_atoms => |payload| {
            if (payload.pivot.index() >= coordinates.len) return error.InvalidAtomIndex;
            const pivot = coordinates[payload.pivot.index()];
            for (affected) |atom| coordinates[atom.index()] = add(pivot, scale(subtract(coordinates[atom.index()], pivot), 0.4));
        },
        .scale_fragment => {
            const factor = stateFactor(1.4, dof.state.current);
            for (affected) |atom| coordinates[atom.index()] = scale(coordinates[atom.index()], factor);
        },
        .invert_bond => |payload| {
            if (payload.pivot.index() >= coordinates.len or payload.bound.index() >= coordinates.len) return error.InvalidAtomIndex;
            const pivot = coordinates[payload.pivot.index()];
            const bond = subtract(coordinates[payload.bound.index()], pivot);
            const line_a = add(pivot, .{ .x = bond.y, .y = -bond.x });
            const line_b = subtract(pivot, .{ .x = bond.y, .y = -bond.x });
            for (affected) |atom| coordinates[atom.index()] = reflect(coordinates[atom.index()], line_a, line_b);
        },
        .flip_ring => |payload| {
            if (payload.pivot_a.index() >= coordinates.len or payload.pivot_b.index() >= coordinates.len) return error.InvalidAtomIndex;
            const line_a = coordinates[payload.pivot_a.index()];
            const line_b = coordinates[payload.pivot_b.index()];
            for (affected) |atom| coordinates[atom.index()] = reflect(coordinates[atom.index()], line_a, line_b);
        },
    }
}

pub fn penalty(dof: core.dof.Dof, affected_count: usize, context: PenaltyContext) !f32 {
    try dof.state.validate();
    const state_level = @as(f32, @floatFromInt((dof.state.current + 1) / 2));
    return switch (dof.payload) {
        .flip_fragment => @as(f32, @floatFromInt(@as(u32, @intFromBool(context.chain_with_chain_parent)) * 10 +
            @as(u32, @intFromBool(dof.state.current != 0 and context.constrained_flip)) * 1000)),
        .change_parent_bond_length => if (dof.state.current == 0) 0 else 200 * state_level,
        .rotate_fragment => if (dof.state.current == 0) 0 else 400 * state_level,
        .scale_atoms => if (dof.state.current == 0) 0 else 50 * @as(f32, @floatFromInt(affected_count)),
        .scale_fragment => if (dof.state.current == 0) 0 else 500 * state_level,
        .invert_bond => if (dof.state.current == 0) 0 else 100,
        .flip_ring => |payload| if (dof.state.current == 0) 0 else 200 * @as(f32, @floatFromInt(payload.penalty_multiplier)),
    };
}

fn stateFactor(base: f32, state: u32) f32 {
    var factor = std.math.pow(f32, base, @floatFromInt((state + 1) / 2));
    if (state % 2 == 0) factor = 1 / factor;
    return factor;
}

fn reflect(point: core.math.Vec2, line_a: core.math.Vec2, line_b: core.math.Vec2) core.math.Vec2 {
    const direction = subtract(line_b, line_a);
    const denominator = direction.x * direction.x + direction.y * direction.y;
    if (denominator == 0) return point;
    const relative = subtract(point, line_a);
    const projection = (relative.x * direction.x + relative.y * direction.y) / denominator;
    const on_line = add(line_a, scale(direction, projection));
    return subtract(scale(on_line, 2), point);
}

fn add(left: core.math.Vec2, right: core.math.Vec2) core.math.Vec2 {
    return .{ .x = left.x + right.x, .y = left.y + right.y };
}

fn subtract(left: core.math.Vec2, right: core.math.Vec2) core.math.Vec2 {
    return .{ .x = left.x - right.x, .y = left.y - right.y };
}

fn scale(value: core.math.Vec2, factor: f32) core.math.Vec2 {
    return .{ .x = value.x * factor, .y = value.y * factor };
}

fn testDof(payload: core.dof.Payload, state: u32, count: u32) core.dof.Dof {
    return .{
        .id = core.ids.DofId.fromIndex(0),
        .fragment = core.ids.FragmentId.fromIndex(0),
        .state = .{ .current = state, .optimal = 0, .count = count, .tier = 0 },
        .payload = payload,
    };
}

test "all seven DOF applications preserve pinned local transforms" {
    const affected = [_]core.ids.AtomId{core.ids.AtomId.fromIndex(2)};

    var coordinates = [_]core.math.Vec2{ .{}, .{}, .{ .x = 2, .y = 3 } };
    try apply(testDof(.{ .flip_fragment = .{} }, 1, 2), &affected, &coordinates);
    try std.testing.expectEqual(core.math.Vec2{ .x = 2, .y = -3 }, coordinates[2]);

    coordinates[2] = .{ .x = 2, .y = 3 };
    try apply(testDof(.{ .scale_fragment = .{} }, 1, 5), &affected, &coordinates);
    try std.testing.expectApproxEqAbs(@as(f32, 2.8), coordinates[2].x, 0.0001);

    coordinates[1] = .{ .x = 1, .y = 1 };
    coordinates[2] = .{ .x = 6, .y = 1 };
    try apply(testDof(.{ .scale_atoms = .{ .pivot = core.ids.AtomId.fromIndex(1) } }, 1, 2), &affected, &coordinates);
    try std.testing.expectEqual(core.math.Vec2{ .x = 3, .y = 1 }, coordinates[2]);

    coordinates[2] = .{ .x = 0, .y = 0 };
    try apply(testDof(.{ .change_parent_bond_length = .{} }, 1, 7), &affected, &coordinates);
    try std.testing.expectApproxEqAbs(@as(f32, 30), coordinates[2].x, 0.0001);

    coordinates[2] = .{ .x = 0, .y = 0 };
    try apply(testDof(.{ .rotate_fragment = .{} }, 1, 5), &affected, &coordinates);
    try std.testing.expectApproxEqAbs(@as(f32, -1.7037087), coordinates[2].x, 0.0001);
    try std.testing.expectApproxEqAbs(@as(f32, -12.940952), coordinates[2].y, 0.0001);

    coordinates[0] = .{};
    coordinates[1] = .{ .x = 1 };
    coordinates[2] = .{ .x = 2, .y = 3 };
    try apply(testDof(.{ .invert_bond = .{ .pivot = core.ids.AtomId.fromIndex(0), .bound = core.ids.AtomId.fromIndex(1) } }, 1, 2), &affected, &coordinates);
    try std.testing.expectEqual(core.math.Vec2{ .x = -2, .y = 3 }, coordinates[2]);

    coordinates[0] = .{};
    coordinates[1] = .{ .x = 1 };
    coordinates[2] = .{ .x = 2, .y = 3 };
    try apply(testDof(.{ .flip_ring = .{
        .ring = core.ids.RingId.fromIndex(0),
        .pivot_a = core.ids.AtomId.fromIndex(0),
        .pivot_b = core.ids.AtomId.fromIndex(1),
        .penalty_multiplier = 3,
    } }, 1, 2), &affected, &coordinates);
    try std.testing.expectEqual(core.math.Vec2{ .x = 2, .y = -3 }, coordinates[2]);
}

test "DOF penalties preserve state levels and flip context" {
    try std.testing.expectEqual(@as(f32, 1010), try penalty(testDof(.{ .flip_fragment = .{} }, 1, 2), 0, .{
        .constrained_flip = true,
        .chain_with_chain_parent = true,
    }));
    try std.testing.expectEqual(@as(f32, 400), try penalty(testDof(.{ .change_parent_bond_length = .{} }, 3, 7), 0, .{}));
    try std.testing.expectEqual(@as(f32, 800), try penalty(testDof(.{ .rotate_fragment = .{} }, 3, 5), 0, .{}));
    try std.testing.expectEqual(@as(f32, 150), try penalty(testDof(.{ .scale_atoms = .{ .pivot = core.ids.AtomId.fromIndex(0) } }, 1, 2), 3, .{}));
    try std.testing.expectEqual(@as(f32, 1000), try penalty(testDof(.{ .scale_fragment = .{} }, 3, 5), 0, .{}));
    try std.testing.expectEqual(@as(f32, 100), try penalty(testDof(.{ .invert_bond = .{ .pivot = core.ids.AtomId.fromIndex(0), .bound = core.ids.AtomId.fromIndex(1) } }, 1, 2), 1, .{}));
    try std.testing.expectEqual(@as(f32, 600), try penalty(testDof(.{ .flip_ring = .{
        .ring = core.ids.RingId.fromIndex(0),
        .pivot_a = core.ids.AtomId.fromIndex(0),
        .pivot_b = core.ids.AtomId.fromIndex(1),
        .penalty_multiplier = 3,
    } }, 1, 2), 0, .{}));
}

const ScoreContext = struct {
    calls: usize = 0,
};

fn targetScore(raw_context: ?*anyopaque, states: []const u32) !f32 {
    const context: *ScoreContext = @ptrCast(@alignCast(raw_context.?));
    context.calls += 1;
    var score_value: f32 = 0;
    for (states) |state| score_value += @as(f32, @floatFromInt(1 -| state)) * 20;
    return score_value;
}

test "solution cache is deterministic, bounded, and chooses lexicographic ties" {
    var context = ScoreContext{};
    var cache = try SolutionCache.init(std.testing.allocator, 2, 1, .{ .context = &context, .scoreFn = targetScore });
    defer cache.deinit();
    try std.testing.expectEqual(@as(f32, 40), try cache.score(&.{ 0, 0 }));
    try std.testing.expectEqual(@as(f32, 40), try cache.score(&.{ 0, 0 }));
    try std.testing.expectEqual(@as(usize, 1), context.calls);
    try std.testing.expectEqual(@as(f32, 20), try cache.score(&.{ 1, 0 }));
    try std.testing.expectEqual(@as(f32, 20), try cache.score(&.{ 0, 1 }));
    const best = try cache.best();
    try std.testing.expectEqualSlices(u32, &.{ 0, 1 }, best.states);

    var bounded_context = ScoreContext{};
    var bounded = try SolutionCache.init(std.testing.allocator, 1, 1, .{ .context = &bounded_context, .scoreFn = targetScore });
    defer bounded.deinit();
    bounded.max_solutions = 0;
    _ = try bounded.score(&.{0});
    try std.testing.expectEqual(rejected_solution_score, try bounded.score(&.{1}));
    try std.testing.expectEqual(@as(usize, 1), bounded_context.calls);
}

test "exhaustive search stores the first clean solution as optimal" {
    var context = ScoreContext{};
    var cache = try SolutionCache.init(std.testing.allocator, 2, 1, .{ .context = &context, .scoreFn = targetScore });
    defer cache.deinit();
    var dofs = [_]core.dof.Dof{
        testDof(.{ .flip_fragment = .{} }, 0, 2),
        testDof(.{ .invert_bond = .{ .pivot = core.ids.AtomId.fromIndex(0), .bound = core.ids.AtomId.fromIndex(1) } }, 0, 2),
    };
    dofs[1].id = core.ids.DofId.fromIndex(1);
    const result = try exhaustiveSearch(std.testing.allocator, &dofs, &.{ 0, 1 }, &cache, 100);
    try std.testing.expect(result.clean_pose);
    try std.testing.expectEqual(@as(f32, 0), result.score);
    try std.testing.expectEqual(@as(u32, 1), dofs[0].state.current);
    try std.testing.expectEqual(@as(u32, 1), dofs[1].state.current);
    try std.testing.expectEqual(@as(usize, 4), context.calls);
}

test "tuple construction retains pinned ascending combination order" {
    var tuples = try buildTuples(std.testing.allocator, 4, 2);
    defer tuples.deinit();
    try std.testing.expectEqual(@as(usize, 6), tuples.count());
    const expected = [_][2]usize{ .{ 0, 1 }, .{ 0, 2 }, .{ 0, 3 }, .{ 1, 2 }, .{ 1, 3 }, .{ 2, 3 } };
    for (expected, 0..) |tuple, index| try std.testing.expectEqualSlices(usize, &tuple, tuples.tuple(index));
}

fn allocateSearchState(allocator: std.mem.Allocator) !void {
    var context = ScoreContext{};
    var cache = try SolutionCache.init(allocator, 2, 1, .{ .context = &context, .scoreFn = targetScore });
    defer cache.deinit();
    var dofs = [_]core.dof.Dof{
        testDof(.{ .flip_fragment = .{} }, 0, 2),
        testDof(.{ .flip_fragment = .{} }, 0, 2),
    };
    _ = try exhaustiveSearch(allocator, &dofs, &.{ 0, 1 }, &cache, 100);
    var tuples = try buildTuples(allocator, 5, 3);
    defer tuples.deinit();
}

test "solution bookkeeping and tuple construction clean every allocation failure" {
    try std.testing.checkAllAllocationFailures(std.testing.allocator, allocateSearchState, .{});
}
