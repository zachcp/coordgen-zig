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
pub const standard_crossing_bond_penalty: f32 = 2500;
pub const terminal_bond_crossing_multiplier: f32 = 0.5;
pub const macrocycle_bond_crossing_multiplier: f32 = 8;
pub const ring_bond_crossing_multiplier: f32 = 2;

pub const Evaluator = struct {
    context: ?*anyopaque = null,
    scoreFn: *const fn (?*anyopaque, []const u32) core.errors.Error!f32,

    fn score(self: Evaluator, states: []const u32) !f32 {
        return self.scoreFn(self.context, states);
    }
};

pub const PoseEvaluator = struct {
    dofs: []core.dof.Dof,
    affected_atoms: []const core.ids.AtomId,
    base_local_coordinates: []const core.math.Vec2,
    local_coordinates: []core.math.Vec2,
    global_coordinates: []core.math.Vec2,
    rebuild_context: ?*anyopaque = null,
    rebuildFn: *const fn (?*anyopaque, []const core.math.Vec2, []core.math.Vec2) core.errors.Error!void,
    score_context: ?*anyopaque = null,
    scorePoseFn: *const fn (?*anyopaque, []const core.math.Vec2, []const core.dof.Dof) core.errors.Error!f32,

    pub fn evaluator(self: *PoseEvaluator) Evaluator {
        return .{ .context = self, .scoreFn = evaluateOpaque };
    }

    fn evaluateOpaque(raw_context: ?*anyopaque, states: []const u32) !f32 {
        const self: *PoseEvaluator = @ptrCast(@alignCast(raw_context orelse return error.InvalidMapping));
        return self.evaluate(states);
    }

    pub fn evaluate(self: *PoseEvaluator, states: []const u32) core.errors.Error!f32 {
        if (states.len != self.dofs.len or
            self.base_local_coordinates.len != self.local_coordinates.len or
            self.local_coordinates.len != self.global_coordinates.len)
        {
            return error.InvalidMapping;
        }
        @memcpy(self.local_coordinates, self.base_local_coordinates);
        for (self.dofs, states) |*dof, state| {
            if (state >= dof.state.count) return error.InvalidMapping;
            dof.state.current = state;
            const affected = try affectedAtoms(self.affected_atoms, dof.*);
            apply(dof.*, affected, self.local_coordinates) catch |err| switch (err) {
                error.InvalidDofState => return error.InvalidMapping,
                else => return err,
            };
        }
        try self.rebuildFn(self.rebuild_context, self.local_coordinates, self.global_coordinates);
        return self.scorePoseFn(self.score_context, self.global_coordinates, self.dofs);
    }
};

pub const FramePose = struct {
    frames: []const core.dof.FragmentFrame,
    rebuild_order: []const core.ids.FragmentId,
    fragment_atoms: []const core.ids.AtomId,
    child_attachments: []const core.dof.ChildAttachment,
    atom_coordinates: []core.math.Vec2,
    attachment_coordinates: []core.math.Vec2,
    global_coordinates: []core.math.Vec2,

    pub fn validate(self: FramePose) core.errors.Error!void {
        if (self.atom_coordinates.len != self.fragment_atoms.len or
            self.attachment_coordinates.len != self.child_attachments.len or
            self.rebuild_order.len != self.frames.len)
        {
            return error.InvalidMapping;
        }
        for (self.frames, 0..) |frame, index| {
            if (frame.id.index() != index or
                rangeEnd(frame.atoms) > self.fragment_atoms.len or
                rangeEnd(frame.attachments) > self.child_attachments.len)
            {
                return error.InvalidMapping;
            }
        }
    }

    pub fn rebuild(self: FramePose) core.errors.Error!void {
        try self.validate();
        for (self.rebuild_order) |fragment_id| {
            if (!fragment_id.isValid() or fragment_id.index() >= self.frames.len) return error.InvalidMapping;
            const frame = self.frames[fragment_id.index()];
            var position: core.math.Vec2 = .{};
            var angle: f32 = 0;
            if (frame.parent.isValid()) {
                const parent_position = try scoreCoordinate(self.global_coordinates, frame.parent_atom);
                position = try scoreCoordinate(self.global_coordinates, frame.anchor_atom);
                const direction = subtract(position, parent_position);
                angle = std.math.atan2(-direction.y, direction.x);
            }
            const sine = @sin(angle);
            const cosine = @cos(angle);
            for (self.fragment_atoms[frame.atoms.start..][0..frame.atoms.len], self.atom_coordinates[frame.atoms.start..][0..frame.atoms.len]) |atom, local| {
                if (!atom.isValid() or atom.index() >= self.global_coordinates.len) return error.InvalidAtomIndex;
                self.global_coordinates[atom.index()] = add(geometry.rotate(local, sine, cosine), position);
            }
            for (self.child_attachments[frame.attachments.start..][0..frame.attachments.len], self.attachment_coordinates[frame.attachments.start..][0..frame.attachments.len]) |attachment, local| {
                if (!attachment.atom.isValid() or attachment.atom.index() >= self.global_coordinates.len) return error.InvalidAtomIndex;
                self.global_coordinates[attachment.atom.index()] = add(geometry.rotate(local, sine, cosine), position);
            }
        }
    }

    fn coordinate(self: FramePose, fragment: core.ids.FragmentId, atom: core.ids.AtomId) core.errors.Error!*core.math.Vec2 {
        if (!fragment.isValid() or fragment.index() >= self.frames.len or !atom.isValid()) return error.InvalidMapping;
        const frame = self.frames[fragment.index()];
        for (self.fragment_atoms[frame.atoms.start..][0..frame.atoms.len], frame.atoms.start..) |candidate, index| {
            if (candidate == atom) return &self.atom_coordinates[index];
        }
        for (self.child_attachments[frame.attachments.start..][0..frame.attachments.len], frame.attachments.start..) |attachment, index| {
            if (attachment.atom == atom) return &self.attachment_coordinates[index];
        }
        return error.InvalidMapping;
    }
};

fn rangeEnd(range: core.dof.AtomRange) usize {
    return @as(usize, range.start) + @as(usize, range.len);
}

pub const FramePoseEvaluator = struct {
    dofs: []core.dof.Dof,
    affected_atoms: []const core.ids.AtomId,
    pose: FramePose,
    base_atom_coordinates: []const core.math.Vec2,
    base_attachment_coordinates: []const core.math.Vec2,
    score_context: ?*anyopaque = null,
    scorePoseFn: *const fn (?*anyopaque, []const core.math.Vec2, []const core.dof.Dof) core.errors.Error!f32,

    pub fn evaluator(self: *FramePoseEvaluator) Evaluator {
        return .{ .context = self, .scoreFn = evaluateOpaque };
    }

    fn evaluateOpaque(raw_context: ?*anyopaque, states: []const u32) !f32 {
        const self: *FramePoseEvaluator = @ptrCast(@alignCast(raw_context orelse return error.InvalidMapping));
        return self.evaluate(states);
    }

    pub fn evaluate(self: *FramePoseEvaluator, states: []const u32) core.errors.Error!f32 {
        if (states.len != self.dofs.len or
            self.base_atom_coordinates.len != self.pose.atom_coordinates.len or
            self.base_attachment_coordinates.len != self.pose.attachment_coordinates.len)
        {
            return error.InvalidMapping;
        }
        @memcpy(self.pose.atom_coordinates, self.base_atom_coordinates);
        @memcpy(self.pose.attachment_coordinates, self.base_attachment_coordinates);
        for (self.dofs, states) |*dof, state| {
            if (state >= dof.state.count) return error.InvalidMapping;
            dof.state.current = state;
            try applyToFrame(dof.*, try affectedAtoms(self.affected_atoms, dof.*), self.pose);
        }
        try self.pose.rebuild();
        // Fragment construction can propagate a scale-atoms DOF onto the
        // parent-side atom of its attachment bond. Upstream stores that atom
        // pointer outside the fragment's local-coordinate map, so the DOF
        // mutates its already-built global coordinate. Apply those external
        // entries after rebuilding; child-side attachments are present in
        // the owner's frame and were handled above.
        for (self.pose.rebuild_order) |fragment_id| {
            for (self.dofs) |dof| {
                if (dof.fragment != fragment_id or dof.state.current == 0) continue;
                const payload = switch (dof.payload) {
                    .scale_atoms => |payload| payload,
                    else => continue,
                };
                const pivot = (try self.pose.coordinate(dof.fragment, payload.pivot)).*;
                for (try affectedAtoms(self.affected_atoms, dof)) |atom| {
                    _ = self.pose.coordinate(dof.fragment, atom) catch |err| switch (err) {
                        error.InvalidMapping => {
                            if (!atom.isValid() or atom.index() >= self.pose.global_coordinates.len) return error.InvalidAtomIndex;
                            const coordinate = &self.pose.global_coordinates[atom.index()];
                            coordinate.* = add(pivot, scale(subtract(coordinate.*, pivot), 0.4));
                            continue;
                        },
                        else => return err,
                    };
                }
            }
        }
        return self.scorePoseFn(self.score_context, self.pose.global_coordinates, self.dofs);
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

    pub fn has(self: SolutionCache, candidate: []const u32) bool {
        if (candidate.len != self.state_count) return false;
        for (0..self.count()) |index| {
            if (std.mem.eql(u32, self.solution(index), candidate)) return true;
        }
        return false;
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

pub const SearchResult = struct {
    score: f32,
    clean_pose: bool,
    iterations: usize,
};

pub fn tieredSearch(
    allocator: std.mem.Allocator,
    dofs: []core.dof.Dof,
    cache: *SolutionCache,
    precision: f32,
    initial_tier: u32,
) core.errors.Error!SearchResult {
    if (cache.state_count != dofs.len or !std.math.isFinite(precision) or precision <= 0) return error.InvalidOption;
    for (dofs) |dof| dof.state.validate() catch return error.InvalidMapping;
    const scaled_width = @as(f64, 6) * @as(f64, precision);
    if (scaled_width > @as(f64, @floatFromInt(std.math.maxInt(usize)))) return error.TooManyItems;
    const max_width = @max(@as(usize, 1), @as(usize, @intFromFloat(scaled_width)));
    var current_tier = initial_tier;
    var growing = CandidateBuffer.init(allocator, dofs.len);
    defer growing.deinit();
    const states = allocator.alloc(u32, dofs.len) catch return error.OutOfMemory;
    defer allocator.free(states);
    saveStates(dofs, states);
    var best_score = try cache.score(states);
    try growing.put(states, best_score);

    var iterations: usize = 0;
    var valid = true;
    while (valid and growing.count() != 0 and iterations < 100) {
        iterations += 1;
        valid = try growSolutions(allocator, dofs, cache, &current_tier, &growing, &best_score, max_width);
    }
    const best = try cache.best();
    loadStates(dofs, best.states);
    // Materialize the selected pose; the final scored candidate is not
    // necessarily the best cached solution.
    _ = try cache.evaluator.score(best.states);
    return .{ .score = best.score, .clean_pose = best.score < clash_energy_threshold, .iterations = iterations };
}

const CandidateBuffer = struct {
    allocator: std.mem.Allocator,
    state_count: usize,
    states: std.ArrayList(u32) = .empty,
    scores: std.ArrayList(f32) = .empty,

    fn init(allocator: std.mem.Allocator, state_count: usize) CandidateBuffer {
        return .{ .allocator = allocator, .state_count = state_count };
    }

    fn deinit(self: *CandidateBuffer) void {
        self.scores.deinit(self.allocator);
        self.states.deinit(self.allocator);
    }

    fn count(self: CandidateBuffer) usize {
        return self.scores.items.len;
    }

    fn candidate(self: CandidateBuffer, index: usize) []const u32 {
        return self.states.items[index * self.state_count ..][0..self.state_count];
    }

    fn clear(self: *CandidateBuffer) void {
        self.states.clearRetainingCapacity();
        self.scores.clearRetainingCapacity();
    }

    fn put(self: *CandidateBuffer, states: []const u32, score_value: f32) !void {
        if (states.len != self.state_count) return error.InvalidMapping;
        for (0..self.count()) |index| {
            if (std.mem.eql(u32, self.candidate(index), states)) {
                self.scores.items[index] = score_value;
                return;
            }
        }
        try self.states.ensureUnusedCapacity(self.allocator, states.len);
        try self.scores.ensureUnusedCapacity(self.allocator, 1);
        self.states.appendSliceAssumeCapacity(states);
        self.scores.appendAssumeCapacity(score_value);
    }

    fn cloneFrom(self: *CandidateBuffer, source: CandidateBuffer) !void {
        self.clear();
        try self.states.appendSlice(self.allocator, source.states.items);
        errdefer self.states.clearRetainingCapacity();
        try self.scores.appendSlice(self.allocator, source.scores.items);
    }
};

fn growSolutions(
    allocator: std.mem.Allocator,
    dofs: []core.dof.Dof,
    cache: *SolutionCache,
    current_tier: *u32,
    growing: *CandidateBuffer,
    best_score: *f32,
    max_width: usize,
) !bool {
    var previous = CandidateBuffer.init(allocator, dofs.len);
    defer previous.deinit();
    try previous.cloneFrom(growing.*);
    const order = allocator.alloc(usize, previous.count()) catch return error.OutOfMemory;
    defer allocator.free(order);
    for (order, 0..) |*index, value| index.* = value;
    sortCandidates(order, previous);
    growing.clear();
    const best_score_for_run = best_score.*;
    const states = allocator.alloc(u32, dofs.len) catch return error.OutOfMemory;
    defer allocator.free(states);

    for (order, 0..) |candidate_index, rank| {
        if (rank > max_width) break;
        for (dofs, 0..) |*dof, dof_index| {
            if (dof.state.tier > current_tier.*) continue;
            loadStates(dofs, previous.candidate(candidate_index));
            for (1..dof.state.count) |_| {
                dofs[dof_index].state.advance();
                saveStates(dofs, states);
                if (cache.has(states)) continue;
                const score_value = try cache.score(states);
                if (score_value == rejected_solution_score) return false;
                if (score_value < best_score.*) best_score.* = score_value;
                if (score_value < best_score_for_run) try growing.put(states, score_value);
            }
        }
    }
    if (growing.count() == 0 and current_tier.* < core.dof.Tier.scale_fragment) {
        current_tier.* += 3;
        try growing.cloneFrom(previous);
    }
    return true;
}

fn sortCandidates(order: []usize, candidates: CandidateBuffer) void {
    for (1..order.len) |index| {
        const value = order[index];
        var cursor = index;
        while (cursor > 0 and candidateLess(candidates, value, order[cursor - 1])) : (cursor -= 1) {
            order[cursor] = order[cursor - 1];
        }
        order[cursor] = value;
    }
}

fn candidateLess(candidates: CandidateBuffer, left: usize, right: usize) bool {
    const left_score = candidates.scores.items[left];
    const right_score = candidates.scores.items[right];
    return left_score < right_score or
        (left_score == right_score and lexicographicLess(candidates.candidate(left), candidates.candidate(right)));
}

fn saveStates(dofs: []const core.dof.Dof, states: []u32) void {
    std.debug.assert(dofs.len == states.len);
    for (dofs, states) |dof, *state| state.* = dof.state.current;
}

fn loadStates(dofs: []core.dof.Dof, states: []const u32) void {
    std.debug.assert(dofs.len == states.len);
    for (dofs, states) |*dof, state| dof.state.current = state;
}

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

pub const BondScoreView = struct {
    start: core.ids.AtomId,
    end: core.ids.AtomId,
    component: core.ids.MoleculeId = .invalid,
    crossing_penalty_multiplier: f32 = 1,
    terminal: bool = false,
    macrocycle: bool = false,
    small_ring: bool = false,
    residue_interaction: bool = false,
};

pub const ResidueInteractionScoreView = struct {
    residue: core.ids.AtomId,
    end: core.ids.AtomId,
    addition_vector: core.math.Vec2 = .{},
};

pub const RingScoreView = struct {
    atoms: []const core.ids.AtomId,
    fragment: core.ids.FragmentId,
    component: core.ids.MoleculeId = .invalid,
};

pub const ProximityScoreView = struct {
    local_molecule: core.ids.MoleculeId,
    local_fragment: core.ids.FragmentId,
    other_molecule: core.ids.MoleculeId,
    addition_vector: core.math.Vec2,
};

pub fn scoreClashInteractions(
    interactions: []const core.interaction.Interaction,
    coordinates: []const core.math.Vec2,
) core.errors.Error!f32 {
    var energy: f32 = 0;
    for (interactions) |interaction| switch (interaction.payload) {
        .clash => |clash| {
            const start = try scoreCoordinate(coordinates, clash.segment_start);
            const point = try scoreCoordinate(coordinates, clash.point);
            const end = try scoreCoordinate(coordinates, clash.segment_end);
            const distance = geometry.squaredDistancePointSegment(point, start, end).squared_distance;
            if (distance <= clash.rest_squared_distance) {
                energy += 0.5 * clash.force_constant * clash.secondary_force_constant * (clash.rest_squared_distance - distance);
            }
        },
        else => {},
    };
    return energy;
}

pub fn scoreCrossBonds(
    bonds: []const BondScoreView,
    residue_bonds: []const BondScoreView,
    residue_interactions: []const ResidueInteractionScoreView,
    score_residue_interactions: bool,
    coordinates: []const core.math.Vec2,
) core.errors.Error!f32 {
    var energy: f32 = 0;
    if (bonds.len > 2) {
        for (bonds, 0..) |first, first_index| {
            if (first.residue_interaction) continue;
            try validateScoreBond(first, coordinates.len);
            for (bonds[first_index + 1 ..]) |second| {
                if (second.residue_interaction) continue;
                try validateScoreBond(second, coordinates.len);
                if (first.component.isValid() and second.component.isValid() and first.component != second.component) continue;
                if (!bondsClash(first, second, coordinates)) continue;
                var penalty_value = standard_crossing_bond_penalty * first.crossing_penalty_multiplier * second.crossing_penalty_multiplier;
                if (first.terminal or second.terminal) penalty_value *= terminal_bond_crossing_multiplier;
                if (first.macrocycle or second.macrocycle) penalty_value *= macrocycle_bond_crossing_multiplier;
                if (first.small_ring or second.small_ring) penalty_value *= ring_bond_crossing_multiplier;
                energy += penalty_value;
            }
        }
    }
    if (score_residue_interactions) {
        energy += try scoreResidueInteractionCrossings(residue_bonds, residue_interactions, coordinates);
    }
    return energy;
}

fn scoreResidueInteractionCrossings(
    bonds: []const BondScoreView,
    interactions: []const ResidueInteractionScoreView,
    coordinates: []const core.math.Vec2,
) core.errors.Error!f32 {
    var energy: f32 = 0;
    var group_start: usize = 0;
    while (group_start < interactions.len) {
        var group_end = group_start + 1;
        while (group_end < interactions.len and interactions[group_end].residue == interactions[group_start].residue) : (group_end += 1) {}
        const group = interactions[group_start..group_end];
        for (group) |interaction| {
            if (!interaction.residue.isValid() or !interaction.end.isValid() or interaction.residue.index() >= coordinates.len or
                interaction.end.index() >= coordinates.len or !interaction.addition_vector.isFinite()) return error.InvalidMapping;
        }
        // Preserve upstream's asymmetric indices exactly: ri1 starts at zero
        // and stops before the last item, while ri2 starts at one for every
        // ri1 rather than at ri1 + 1 (CoordgenMinimizer.cpp:871-875).
        if (group.len > 1) for (group[0 .. group.len - 1]) |first| for (group[1..]) |second| {
            const first_end = coordinates[first.end.index()];
            const second_end = coordinates[second.end.index()];
            if (geometry.segmentIntersection(
                geometry.add(first_end, geometry.scale(first.addition_vector, 0.2)),
                geometry.add(second_end, geometry.scale(second.addition_vector, 0.2)),
                first_end,
                second_end,
            ) != null) energy += 15;
            for (bonds) |bond| {
                if (!bond.start.isValid() or !bond.end.isValid() or bond.start.index() >= coordinates.len or
                    bond.end.index() >= coordinates.len) return error.InvalidMapping;
                if (bond.start == first.end or bond.end == first.end or bond.start == second.end or bond.end == second.end) continue;
                if (geometry.segmentIntersection(
                    first_end,
                    second_end,
                    coordinates[bond.start.index()],
                    coordinates[bond.end.index()],
                ) != null) energy += 10;
            }
        };
        group_start = group_end;
    }
    return energy;
}

pub fn scoreAtomsInsideRings(
    rings: []const RingScoreView,
    atom_fragments: []const core.ids.FragmentId,
    fragment_components: []const core.ids.MoleculeId,
    coordinates: []const core.math.Vec2,
) core.errors.Error!f32 {
    if (atom_fragments.len != coordinates.len) return error.InvalidMapping;
    var energy: f32 = 0;
    for (rings) |ring| {
        if (ring.atoms.len < 3 or ring.atoms.len >= 9) continue;
        var center: core.math.Vec2 = .{};
        for (ring.atoms) |atom| center = add(center, try scoreCoordinate(coordinates, atom));
        center = scale(center, 1 / @as(f32, @floatFromInt(ring.atoms.len)));
        for (coordinates, atom_fragments) |coordinate_value, fragment| {
            if (fragment == ring.fragment) continue;
            if (fragment_components.len != 0) {
                if (!fragment.isValid() or fragment.index() >= fragment_components.len) return error.InvalidMapping;
                if (ring.component.isValid() and fragment_components[fragment.index()] != ring.component) continue;
            }
            const difference = subtract(center, coordinate_value);
            if (difference.x > core.math.bond_length or difference.y > core.math.bond_length or
                difference.x < -core.math.bond_length or difference.y < -core.math.bond_length)
            {
                continue;
            }
            const squared = difference.x * difference.x + difference.y * difference.y;
            if (squared > core.math.bond_length * core.math.bond_length) continue;
            const distance = @sqrt(squared);
            if (distance < core.math.bond_length) energy += 50 + 100 * (1 - distance / core.math.bond_length);
        }
    }
    return energy;
}

pub fn scoreProximityRelationsOnOppositeSides(relations: []const ProximityScoreView) core.errors.Error!f32 {
    var energy: f32 = 0;
    for (relations, 0..) |first, first_index| {
        try validateProximity(first);
        if (first.other_molecule == first.local_molecule) continue;
        for (relations[first_index + 1 ..]) |second| {
            try validateProximity(second);
            if (second.local_molecule != first.local_molecule or
                second.other_molecule == second.local_molecule or
                second.local_fragment == first.local_fragment or
                second.other_molecule != first.other_molecule)
            {
                continue;
            }
            const angle = geometry.unsignedAngle(first.addition_vector, .{}, second.addition_vector);
            if (angle > 90) energy += 100 + 50 * (angle - 90);
        }
    }
    return energy;
}

fn validateProximity(relation: ProximityScoreView) core.errors.Error!void {
    if (!relation.local_molecule.isValid() or !relation.local_fragment.isValid() or
        !relation.other_molecule.isValid() or !relation.addition_vector.isFinite())
    {
        return error.InvalidMapping;
    }
}

/// Apply upstream's terminal-bond fallback after discrete search fails. The
/// caller re-scores the pose after this deterministic coordinate mutation.
pub fn avoidTerminalClashes(
    initial_score: f32,
    bonds: []const BondScoreView,
    atom_degrees: []const u32,
    fixed: []const bool,
    coordinates: []core.math.Vec2,
) core.errors.Error!bool {
    if (atom_degrees.len != coordinates.len or fixed.len != coordinates.len) return error.InvalidMapping;
    if (initial_score < 0.1) return false;
    var changed = false;
    for (bonds) |bond| {
        if (bond.residue_interaction or !bond.terminal) continue;
        try validateScoreBond(bond, coordinates.len);
        var terminal = bond.end;
        var root = bond.start;
        if (atom_degrees[terminal.index()] != 1) {
            terminal = bond.start;
            root = bond.end;
        }
        if (fixed[terminal.index()]) continue;
        for (bonds) |other| {
            if (other.residue_interaction) continue;
            try validateScoreBond(other, coordinates.len);
            if (!bondsClash(bond, other, coordinates)) continue;
            coordinates[terminal.index()] = add(
                coordinates[root.index()],
                scale(subtract(coordinates[terminal.index()], coordinates[root.index()]), 0.1),
            );
            changed = true;
        }
    }
    return changed;
}

fn validateScoreBond(bond: BondScoreView, atom_count: usize) core.errors.Error!void {
    if (!bond.start.isValid() or !bond.end.isValid() or bond.start.index() >= atom_count or
        bond.end.index() >= atom_count or bond.start == bond.end or
        !std.math.isFinite(bond.crossing_penalty_multiplier))
    {
        return error.InvalidMapping;
    }
}

fn bondsClash(first: BondScoreView, second: BondScoreView, coordinates: []const core.math.Vec2) bool {
    if (first.start == second.start or first.start == second.end or first.end == second.start or first.end == second.end) return false;
    const first_start = coordinates[first.start.index()];
    const first_end = coordinates[first.end.index()];
    const second_start = coordinates[second.start.index()];
    const second_end = coordinates[second.end.index()];
    const coincidence_limit = sketcher_epsilon * sketcher_epsilon;
    if (geometry.squaredDistance(first_start, second_start) < coincidence_limit or
        geometry.squaredDistance(first_start, second_end) < coincidence_limit or
        geometry.squaredDistance(first_end, second_start) < coincidence_limit or
        geometry.squaredDistance(first_end, second_end) < coincidence_limit)
    {
        return true;
    }
    return geometry.segmentIntersection(first_start, first_end, second_start, second_end) != null;
}

fn scoreCoordinate(coordinates: []const core.math.Vec2, atom: core.ids.AtomId) core.errors.Error!core.math.Vec2 {
    if (!atom.isValid() or atom.index() >= coordinates.len) return error.InvalidAtomIndex;
    return coordinates[atom.index()];
}

/// The atoms one DOF moves, as a checked slice of the collection's flat array.
///
/// `core.dof.Dof` stores a start/len range rather than a slice so no DOF owns
/// pointers, which means every consumer has to bounds-check the range against
/// the array it came from. This is that check, in one place: it used to be
/// open-coded identically at both evaluator sites while this function had no
/// callers at all, which is how the two copies could have drifted apart
/// unnoticed (cgz-7v2.24).
pub fn affectedAtoms(
    collection_atoms: []const core.ids.AtomId,
    dof: core.dof.Dof,
) core.errors.Error![]const core.ids.AtomId {
    const start = dof.affected_atoms.start;
    const end = std.math.add(u32, start, dof.affected_atoms.len) catch return error.InvalidMapping;
    if (end > collection_atoms.len) return error.InvalidMapping;
    return collection_atoms[start..end];
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

pub fn applyToFrame(dof: core.dof.Dof, affected: []const core.ids.AtomId, pose: FramePose) core.errors.Error!void {
    dof.state.validate() catch return error.InvalidMapping;
    try pose.validate();
    if (dof.state.current == 0) return;
    switch (dof.payload) {
        .flip_fragment => try transformWholeFrame(pose, dof.fragment, .flip_y),
        .change_parent_bond_length => {
            const factor = stateFactor(1.6, dof.state.current);
            try transformWholeFrame(pose, dof.fragment, .{ .translate_x = core.math.bond_length * (factor - 1) });
        },
        .rotate_fragment => {
            var angle = @as(f32, std.math.pi) / 180 * 15 * @as(f32, @floatFromInt((dof.state.current + 1) / 2));
            if (dof.state.current % 2 == 0) angle = -angle;
            try transformWholeFrame(pose, dof.fragment, .{ .rotate = .{
                .origin = .{ .x = -core.math.bond_length },
                .sine = @sin(angle),
                .cosine = @cos(angle),
            } });
        },
        .scale_atoms => |payload| {
            const pivot = (try pose.coordinate(dof.fragment, payload.pivot)).*;
            for (affected) |atom| {
                const coordinate_value = pose.coordinate(dof.fragment, atom) catch |err| switch (err) {
                    // The parent-side attachment atom is not part of this
                    // fragment's local-coordinate map. FramePoseEvaluator
                    // applies that upstream pointer side effect after rebuild.
                    error.InvalidMapping => continue,
                    else => return err,
                };
                coordinate_value.* = add(pivot, scale(subtract(coordinate_value.*, pivot), 0.4));
            }
        },
        .scale_fragment => try transformWholeFrame(pose, dof.fragment, .{ .scale = stateFactor(1.4, dof.state.current) }),
        .invert_bond => |payload| {
            const pivot = (try pose.coordinate(dof.fragment, payload.pivot)).*;
            const bond = subtract((try pose.coordinate(dof.fragment, payload.bound)).*, pivot);
            const line_a = add(pivot, .{ .x = bond.y, .y = -bond.x });
            const line_b = subtract(pivot, .{ .x = bond.y, .y = -bond.x });
            for (affected) |atom| {
                const coordinate_value = try pose.coordinate(dof.fragment, atom);
                coordinate_value.* = reflect(coordinate_value.*, line_a, line_b);
            }
        },
        .flip_ring => |payload| {
            const line_a = (try pose.coordinate(dof.fragment, payload.pivot_a)).*;
            const line_b = (try pose.coordinate(dof.fragment, payload.pivot_b)).*;
            for (affected) |atom| {
                const coordinate_value = try pose.coordinate(dof.fragment, atom);
                coordinate_value.* = reflect(coordinate_value.*, line_a, line_b);
            }
        },
    }
}

const FrameTransform = union(enum) {
    flip_y,
    translate_x: f32,
    rotate: struct { origin: core.math.Vec2, sine: f32, cosine: f32 },
    scale: f32,
};

fn transformWholeFrame(pose: FramePose, fragment: core.ids.FragmentId, transform: FrameTransform) core.errors.Error!void {
    if (!fragment.isValid() or fragment.index() >= pose.frames.len) return error.InvalidMapping;
    const frame = pose.frames[fragment.index()];
    for (pose.atom_coordinates[frame.atoms.start..][0..frame.atoms.len]) |*coordinate_value| applyFrameTransform(coordinate_value, transform);
    for (pose.attachment_coordinates[frame.attachments.start..][0..frame.attachments.len]) |*coordinate_value| applyFrameTransform(coordinate_value, transform);
}

fn applyFrameTransform(coordinate_value: *core.math.Vec2, transform: FrameTransform) void {
    switch (transform) {
        .flip_y => coordinate_value.y = -coordinate_value.y,
        .translate_x => |amount| coordinate_value.x += amount,
        .rotate => |rotation| coordinate_value.* = add(
            geometry.rotate(subtract(coordinate_value.*, rotation.origin), rotation.sine, rotation.cosine),
            rotation.origin,
        ),
        .scale => |factor| coordinate_value.* = scale(coordinate_value.*, factor),
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

test "frame rebuild preserves duplicate child anchors and parent DOFs move the child" {
    const frames = [_]core.dof.FragmentFrame{
        .{
            .id = core.ids.FragmentId.fromIndex(0),
            .atoms = .{ .start = 0, .len = 1 },
            .attachments = .{ .start = 0, .len = 1 },
        },
        .{
            .id = core.ids.FragmentId.fromIndex(1),
            .parent = core.ids.FragmentId.fromIndex(0),
            .parent_atom = core.ids.AtomId.fromIndex(0),
            .anchor_atom = core.ids.AtomId.fromIndex(1),
            .atoms = .{ .start = 1, .len = 2 },
            .attachments = .{ .start = 1, .len = 0 },
        },
    };
    const order = [_]core.ids.FragmentId{ core.ids.FragmentId.fromIndex(0), core.ids.FragmentId.fromIndex(1) };
    const fragment_atoms = [_]core.ids.AtomId{
        core.ids.AtomId.fromIndex(0),
        core.ids.AtomId.fromIndex(1),
        core.ids.AtomId.fromIndex(2),
    };
    const attachments = [_]core.dof.ChildAttachment{.{
        .child = core.ids.FragmentId.fromIndex(1),
        .atom = core.ids.AtomId.fromIndex(1),
    }};
    var local_atoms = [_]core.math.Vec2{ .{}, .{}, .{ .x = 50 } };
    var local_attachments = [_]core.math.Vec2{.{ .x = 50, .y = 20 }};
    var global: [3]core.math.Vec2 = undefined;
    const pose = FramePose{
        .frames = &frames,
        .rebuild_order = &order,
        .fragment_atoms = &fragment_atoms,
        .child_attachments = &attachments,
        .atom_coordinates = &local_atoms,
        .attachment_coordinates = &local_attachments,
        .global_coordinates = &global,
    };
    try pose.rebuild();
    try std.testing.expectApproxEqAbs(@as(f32, 50), global[1].x, 0.0001);
    try std.testing.expectApproxEqAbs(@as(f32, 20), global[1].y, 0.0001);

    try applyToFrame(testDof(.{ .flip_fragment = .{} }, 1, 2), &.{}, pose);
    try pose.rebuild();
    try std.testing.expectEqual(core.math.Vec2{}, local_atoms[1]);
    try std.testing.expectApproxEqAbs(@as(f32, 50), global[1].x, 0.0001);
    try std.testing.expectApproxEqAbs(@as(f32, -20), global[1].y, 0.0001);
    try std.testing.expectApproxEqAbs(@as(f32, 96.423836), global[2].x, 0.0001);
    try std.testing.expectApproxEqAbs(@as(f32, -38.569534), global[2].y, 0.0001);
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

test "discrete clash crossing and atoms-in-rings scores preserve pinned penalties" {
    const coordinates = [_]core.math.Vec2{
        .{ .x = -2 },
        .{ .x = 2 },
        .{ .y = -2 },
        .{ .y = 2 },
        .{ .x = 10 },
        .{ .x = 12 },
    };
    const clash = core.interaction.Interaction{
        .id = core.ids.InteractionId.fromIndex(0),
        .payload = .{ .clash = .{
            .segment_start = core.ids.AtomId.fromIndex(0),
            .point = core.ids.AtomId.fromIndex(2),
            .segment_end = core.ids.AtomId.fromIndex(1),
            .rest_squared_distance = 9,
        } },
    };
    try std.testing.expectApproxEqAbs(@as(f32, 0.25), try scoreClashInteractions(&.{clash}, &coordinates), 0.0001);

    const bonds = [_]BondScoreView{
        .{ .start = core.ids.AtomId.fromIndex(0), .end = core.ids.AtomId.fromIndex(1), .terminal = true, .crossing_penalty_multiplier = 2 },
        .{ .start = core.ids.AtomId.fromIndex(2), .end = core.ids.AtomId.fromIndex(3), .macrocycle = true, .small_ring = true, .crossing_penalty_multiplier = 3 },
        .{ .start = core.ids.AtomId.fromIndex(4), .end = core.ids.AtomId.fromIndex(5) },
    };
    try std.testing.expectEqual(@as(f32, 120_000), try scoreCrossBonds(&bonds, &bonds, &.{}, false, &coordinates));

    const residue_interactions = [_]ResidueInteractionScoreView{
        .{ .residue = core.ids.AtomId.fromIndex(5), .end = core.ids.AtomId.fromIndex(0) },
        .{ .residue = core.ids.AtomId.fromIndex(5), .end = core.ids.AtomId.fromIndex(1) },
    };
    try std.testing.expectEqual(@as(f32, 0), try scoreCrossBonds(&bonds, &bonds, &residue_interactions, false, &coordinates));
    try std.testing.expectEqual(@as(f32, 10), try scoreCrossBonds(&bonds, &bonds, &residue_interactions, true, &coordinates));

    var terminal_coordinates = coordinates;
    const degrees = [_]u32{ 2, 1, 1, 1, 1, 1 };
    const fixed = [_]bool{ false, false, false, false, false, false };
    try std.testing.expect(try avoidTerminalClashes(120_000, &bonds, &degrees, &fixed, &terminal_coordinates));
    try std.testing.expectEqual(core.math.Vec2{ .x = -1.6 }, terminal_coordinates[1]);

    const ring_atoms = [_]core.ids.AtomId{
        core.ids.AtomId.fromIndex(0),
        core.ids.AtomId.fromIndex(1),
        core.ids.AtomId.fromIndex(2),
    };
    const ring_coordinates = [_]core.math.Vec2{ .{ .x = -1 }, .{ .x = 1 }, .{}, .{} };
    const fragments = [_]core.ids.FragmentId{
        core.ids.FragmentId.fromIndex(0),
        core.ids.FragmentId.fromIndex(0),
        core.ids.FragmentId.fromIndex(0),
        core.ids.FragmentId.fromIndex(1),
    };
    try std.testing.expectEqual(@as(f32, 150), try scoreAtomsInsideRings(&.{.{
        .atoms = &ring_atoms,
        .fragment = core.ids.FragmentId.fromIndex(0),
    }}, &fragments, &.{}, &ring_coordinates));

    const proximity = [_]ProximityScoreView{
        .{
            .local_molecule = core.ids.MoleculeId.fromIndex(0),
            .local_fragment = core.ids.FragmentId.fromIndex(0),
            .other_molecule = core.ids.MoleculeId.fromIndex(1),
            .addition_vector = .{ .x = 1 },
        },
        .{
            .local_molecule = core.ids.MoleculeId.fromIndex(0),
            .local_fragment = core.ids.FragmentId.fromIndex(1),
            .other_molecule = core.ids.MoleculeId.fromIndex(1),
            .addition_vector = .{ .x = -1 },
        },
        .{
            .local_molecule = core.ids.MoleculeId.fromIndex(0),
            .local_fragment = core.ids.FragmentId.fromIndex(2),
            .other_molecule = core.ids.MoleculeId.fromIndex(2),
            .addition_vector = .{ .y = 1 },
        },
    };
    try std.testing.expectEqual(@as(f32, 4600), try scoreProximityRelationsOnOppositeSides(&proximity));
}

const RebuildTestContext = struct { calls: usize = 0 };

fn rebuildTestPose(raw_context: ?*anyopaque, local: []const core.math.Vec2, global: []core.math.Vec2) !void {
    const context: *RebuildTestContext = @ptrCast(@alignCast(raw_context.?));
    context.calls += 1;
    for (local, global) |source, *target| target.* = .{ .x = source.x + 10, .y = source.y };
}

fn scoreTestPose(_: ?*anyopaque, coordinates: []const core.math.Vec2, _: []const core.dof.Dof) !f32 {
    return coordinates[2].y;
}

test "pose evaluation resets local coordinates, applies DOFs, and rebuilds once per uncached state" {
    const base = [_]core.math.Vec2{ .{}, .{}, .{ .x = 2, .y = 3 } };
    var local: [3]core.math.Vec2 = undefined;
    var global: [3]core.math.Vec2 = undefined;
    const affected = [_]core.ids.AtomId{core.ids.AtomId.fromIndex(2)};
    var dofs = [_]core.dof.Dof{testDof(.{ .flip_fragment = .{} }, 0, 2)};
    dofs[0].affected_atoms = .{ .start = 0, .len = 1 };
    var rebuild_context = RebuildTestContext{};
    var pose = PoseEvaluator{
        .dofs = &dofs,
        .affected_atoms = &affected,
        .base_local_coordinates = &base,
        .local_coordinates = &local,
        .global_coordinates = &global,
        .rebuild_context = &rebuild_context,
        .rebuildFn = rebuildTestPose,
        .scorePoseFn = scoreTestPose,
    };
    var cache = try SolutionCache.init(std.testing.allocator, 1, 1, pose.evaluator());
    defer cache.deinit();
    try std.testing.expectEqual(@as(f32, -3), try cache.score(&.{1}));
    try std.testing.expectEqual(@as(f32, -3), try cache.score(&.{1}));
    try std.testing.expectEqual(@as(usize, 1), rebuild_context.calls);
    try std.testing.expectEqual(core.math.Vec2{ .x = 12, .y = -3 }, global[2]);
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

test "precision changes the scored-solution budget" {
    var context = ScoreContext{};
    var quick = try SolutionCache.init(std.testing.allocator, 1, 0.2, .{ .context = &context, .scoreFn = targetScore });
    defer quick.deinit();
    var best = try SolutionCache.init(std.testing.allocator, 1, 3, .{ .context = &context, .scoreFn = targetScore });
    defer best.deinit();
    try std.testing.expect(quick.max_solutions < best.max_solutions);
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

test "tiered search widens by three tiers and restores the best solution" {
    var context = ScoreContext{};
    var cache = try SolutionCache.init(std.testing.allocator, 2, 1, .{ .context = &context, .scoreFn = targetScore });
    defer cache.deinit();
    var dofs = [_]core.dof.Dof{
        testDof(.{ .flip_fragment = .{} }, 0, 2),
        testDof(.{ .scale_atoms = .{ .pivot = core.ids.AtomId.fromIndex(0) } }, 0, 2),
    };
    dofs[1].id = core.ids.DofId.fromIndex(1);
    dofs[1].state.tier = core.dof.Tier.scale_atoms;
    const result = try tieredSearch(std.testing.allocator, &dofs, &cache, 1, 0);
    try std.testing.expect(result.clean_pose);
    try std.testing.expectEqual(@as(f32, 0), result.score);
    try std.testing.expectEqual(@as(u32, 1), dofs[0].state.current);
    try std.testing.expectEqual(@as(u32, 1), dofs[1].state.current);
    try std.testing.expect(result.iterations >= 4);
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
    _ = try tieredSearch(allocator, &dofs, &cache, 1, 0);
    var tuples = try buildTuples(allocator, 5, 3);
    defer tuples.deinit();
}

test "solution bookkeeping and tuple construction clean every allocation failure" {
    try core.oom.checkAllocationFailures(std.testing.allocator, allocateSearchState, .{});
}

test "the affected-atom range is bounds-checked rather than trusted" {
    const atoms = [_]core.ids.AtomId{
        core.ids.AtomId.fromIndex(3),
        core.ids.AtomId.fromIndex(4),
    };
    var dof = testDof(.{ .flip_fragment = .{} }, 0, 2);

    dof.affected_atoms = .{ .start = 1, .len = 1 };
    try std.testing.expectEqualSlices(
        core.ids.AtomId,
        atoms[1..2],
        try affectedAtoms(&atoms, dof),
    );

    // Past the end, and an end that overflows u32: both are the collection
    // disagreeing with the range, and neither may produce a slice.
    dof.affected_atoms = .{ .start = 1, .len = 2 };
    try std.testing.expectError(error.InvalidMapping, affectedAtoms(&atoms, dof));
    dof.affected_atoms = .{ .start = std.math.maxInt(u32), .len = 2 };
    try std.testing.expectError(error.InvalidMapping, affectedAtoms(&atoms, dof));
}
