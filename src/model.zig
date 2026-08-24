const std = @import("std");
const core = @import("core");

pub const Atom = struct {
    id: core.ids.AtomId,
    input_index: core.ids.InputIndex,
    atomic_number: core.chemistry.AtomicNumber,
    formal_charge: core.chemistry.FormalCharge = 0,
    fixed: bool = false,
    constrained: bool = false,
    hidden: bool = false,
    template_coordinates: ?core.math.Vec2 = null,
    coordinates_3d: ?core.math.Vec3 = null,
    coordinates: core.math.Vec2 = .{},
    stereo: core.chemistry.AtomStereo = .unspecified,
    stereo_looking_from: core.ids.AtomId = .invalid,
    stereo_atom_a: core.ids.AtomId = .invalid,
    stereo_atom_b: core.ids.AtomId = .invalid,
    cross_layout: bool = false,
};

pub const Bond = struct {
    id: core.ids.BondId,
    input_index: core.ids.InputIndex,
    start: core.ids.AtomId,
    end: core.ids.AtomId,
    input_order: core.chemistry.BondOrder,
    /// Observable normalization result; caller input is never mutated.
    effective_order: core.chemistry.BondOrder,
    skip: bool = false,
    stereo: core.chemistry.BondStereo = .unspecified,
    stereo_atom_a: core.ids.AtomId = .invalid,
    stereo_atom_b: core.ids.AtomId = .invalid,
    display: core.chemistry.BondDisplay = .none,
    crossing_penalty_multiplier: f32 = 1,
};

pub const Ring = struct {
    id: core.ids.RingId,
    atom_start: u32,
    atom_count: u32,
    bond_start: u32,
    bond_count: u32,
};

pub const Fragment = struct {
    id: core.ids.FragmentId,
    parent: core.ids.FragmentId = .invalid,
    atom_start: u32,
    atom_count: u32,
};

pub const Residue = struct {
    id: core.ids.ResidueId,
    atom: core.ids.AtomId,
    chain_start: u32,
    chain_len: u32,
    residue_number: i32,
    closest_ligand_atom: core.ids.AtomId = .invalid,
};

pub const IdRange = extern struct {
    start: u32 = 0,
    len: u32 = 0,
};

pub const ResidueInteraction = struct {
    start: core.ids.AtomId,
    end: core.ids.AtomId,
    other_start_atoms: IdRange = .{},
    other_end_atoms: IdRange = .{},
    crossing_penalty_multiplier: f32 = 1,
};

pub const ExtraBond = struct {
    input_index: core.ids.InputIndex,
    start: core.ids.AtomId,
    end: core.ids.AtomId,
    input_order: core.chemistry.BondOrder = .single,
    effective_order: core.chemistry.BondOrder = .single,
    skip: bool = false,
    stereo: core.chemistry.BondStereo = .unspecified,
    stereo_atom_a: core.ids.AtomId = .invalid,
    stereo_atom_b: core.ids.AtomId = .invalid,
    display: core.chemistry.BondDisplay = .none,
    crossing_penalty_multiplier: f32 = 1,
};

pub const ProximityRelation = union(enum) {
    bond: core.ids.BondId,
    extra_bond: core.ids.InputIndex,
    residue_interaction: core.ids.InputIndex,
};

/// Owns both directions of the canonical/internal permutation. The source
/// slice is `internal_index -> input_index`; construction proves bijection and
/// derives the inverse. No algorithm may infer caller order from storage order.
pub const InputOrderMap = struct {
    allocator: std.mem.Allocator,
    input_to_internal: []u32,
    internal_to_input: []core.ids.InputIndex,

    pub fn init(
        allocator: std.mem.Allocator,
        internal_to_input_source: []const core.ids.InputIndex,
    ) core.errors.Error!InputOrderMap {
        if (internal_to_input_source.len > std.math.maxInt(u32)) return error.TooManyItems;
        const internal_to_input = allocator.dupe(
            core.ids.InputIndex,
            internal_to_input_source,
        ) catch return error.OutOfMemory;
        errdefer allocator.free(internal_to_input);
        const input_to_internal = allocator.alloc(u32, internal_to_input_source.len) catch {
            return error.OutOfMemory;
        };
        errdefer allocator.free(input_to_internal);
        @memset(input_to_internal, std.math.maxInt(u32));

        for (internal_to_input, 0..) |input_index, internal_index| {
            if (input_index >= internal_to_input.len) return error.InvalidMapping;
            if (input_to_internal[input_index] != std.math.maxInt(u32)) {
                return error.InvalidMapping;
            }
            input_to_internal[input_index] = @intCast(internal_index);
        }
        return .{
            .allocator = allocator,
            .input_to_internal = input_to_internal,
            .internal_to_input = internal_to_input,
        };
    }

    pub fn initIdentity(allocator: std.mem.Allocator, count: usize) core.errors.Error!InputOrderMap {
        if (count > std.math.maxInt(u32)) return error.TooManyItems;
        const source = allocator.alloc(core.ids.InputIndex, count) catch return error.OutOfMemory;
        defer allocator.free(source);
        for (source, 0..) |*index, value| index.* = @intCast(value);
        return init(allocator, source);
    }

    pub fn deinit(self: *InputOrderMap) void {
        self.allocator.free(self.internal_to_input);
        self.allocator.free(self.input_to_internal);
        self.* = undefined;
    }

    pub fn writeInputOrder(
        self: InputOrderMap,
        comptime T: type,
        internal_values: []const T,
        output: []T,
    ) core.errors.Error!void {
        if (internal_values.len != self.internal_to_input.len or output.len != internal_values.len) {
            return error.InvalidMapping;
        }
        for (internal_values, self.internal_to_input) |value, input_index| {
            output[input_index] = value;
        }
    }

    /// Upstream removes hidden atoms before exposing canonical order. Internal
    /// storage retains them after `active_count`; these helpers serialize the
    /// exact observable partial permutation without losing owned hidden data.
    pub fn activeInternalToInput(self: InputOrderMap, active_count: u32) []const core.ids.InputIndex {
        std.debug.assert(active_count <= self.internal_to_input.len);
        return self.internal_to_input[0..active_count];
    }

    pub fn writeActiveInputToInternal(
        self: InputOrderMap,
        active_count: u32,
        output: []u32,
    ) core.errors.Error!void {
        if (active_count > self.internal_to_input.len or output.len != self.input_to_internal.len) {
            return error.InvalidMapping;
        }
        @memset(output, std.math.maxInt(u32));
        for (self.internal_to_input[0..active_count], 0..) |input_index, internal_index| {
            output[input_index] = @intCast(internal_index);
        }
    }
};

/// Long-lived native graph storage for one generation context. Phase-local
/// adjacency, rings, fragments, DOFs, and interactions live in their owning
/// arenas/collections rather than hiding allocations inside these values.
pub const WorkingGraph = struct {
    allocator: std.mem.Allocator,
    atoms: []Atom,
    bonds: []Bond,
    extra_bonds: []ExtraBond,
    residues: []Residue,
    residue_interactions: []ResidueInteraction,
    /// Backing storage for both additional endpoint ranges above.
    residue_interaction_atoms: []core.ids.AtomId,
    /// Backing storage for every residue chain range.
    string_bytes: []u8,
    order: InputOrderMap,
    /// Canonical visible atoms occupy this prefix. Hidden atoms remain in the
    /// owned backing array so caller-indexed stereo and bond data stay valid,
    /// but are excluded from every structural phase.
    active_atom_count: u32,
    /// Zero-order ordinary and extra bonds retained outside structural
    /// adjacency in the same stable groups used by upstream initialization.
    proximity_relations: []ProximityRelation,

    pub fn deinit(self: *WorkingGraph) void {
        self.allocator.free(self.proximity_relations);
        self.order.deinit();
        self.allocator.free(self.string_bytes);
        self.allocator.free(self.residue_interaction_atoms);
        self.allocator.free(self.residue_interactions);
        self.allocator.free(self.residues);
        self.allocator.free(self.extra_bonds);
        self.allocator.free(self.bonds);
        self.allocator.free(self.atoms);
        self.* = undefined;
    }
};

test "adversarial internal order emits coordinates in caller input order" {
    const internal_to_input = [_]u32{ 3, 0, 4, 1, 2 };
    var order = try InputOrderMap.init(std.testing.allocator, &internal_to_input);
    defer order.deinit();

    try std.testing.expectEqualSlices(u32, &[_]u32{ 1, 3, 4, 0, 2 }, order.input_to_internal);
    const internal = [_]core.math.Vec2{
        .{ .x = 30, .y = 0 },
        .{ .x = 0, .y = 0 },
        .{ .x = 40, .y = 0 },
        .{ .x = 10, .y = 0 },
        .{ .x = 20, .y = 0 },
    };
    var output: [5]core.math.Vec2 = undefined;
    try order.writeInputOrder(core.math.Vec2, &internal, &output);
    for (output, 0..) |position, input_index| {
        try std.testing.expectEqual(@as(f32, @floatFromInt(input_index * 10)), position.x);
    }
}

test "mapping rejects duplicate and missing input positions" {
    try std.testing.expectError(
        error.InvalidMapping,
        InputOrderMap.init(std.testing.allocator, &[_]u32{ 0, 0 }),
    );
    try std.testing.expectError(
        error.InvalidMapping,
        InputOrderMap.init(std.testing.allocator, &[_]u32{ 0, 2 }),
    );

    var order = try InputOrderMap.initIdentity(std.testing.allocator, 2);
    defer order.deinit();
    var short_output: [1]u32 = undefined;
    try std.testing.expectError(
        error.InvalidMapping,
        order.writeInputOrder(u32, &.{ 10, 20 }, &short_output),
    );
    var full_output: [2]u32 = undefined;
    try std.testing.expectError(
        error.InvalidMapping,
        order.writeInputOrder(u32, &.{10}, &full_output),
    );

    var active_output: [2]u32 = undefined;
    try order.writeActiveInputToInternal(1, &active_output);
    try std.testing.expectEqualSlices(u32, &.{ 0, std.math.maxInt(u32) }, &active_output);
    try std.testing.expectEqualSlices(u32, &.{0}, order.activeInternalToInput(1));
    try std.testing.expectError(
        error.InvalidMapping,
        order.writeActiveInputToInternal(3, &active_output),
    );
}

fn constructAndDiscardOrder(
    allocator: std.mem.Allocator,
    source: []const core.ids.InputIndex,
) !void {
    var order = try InputOrderMap.init(allocator, source);
    order.deinit();
}

fn constructAndDiscardIdentity(allocator: std.mem.Allocator, count: usize) !void {
    var order = try InputOrderMap.initIdentity(allocator, count);
    order.deinit();
}

test "input-order maps report and clean up every allocation failure" {
    const source = [_]core.ids.InputIndex{ 2, 0, 3, 1 };

    var map_counter = std.testing.FailingAllocator.init(std.testing.allocator, .{});
    try constructAndDiscardOrder(map_counter.allocator(), &source);
    // The retained permutation and its inverse are separate allocation sites.
    try std.testing.expectEqual(@as(usize, 2), map_counter.alloc_index);
    try core.oom.checkAllocationFailures(
        std.testing.allocator,
        constructAndDiscardOrder,
        .{&source},
    );

    var identity_counter = std.testing.FailingAllocator.init(std.testing.allocator, .{});
    try constructAndDiscardIdentity(identity_counter.allocator(), source.len);
    // Identity construction also allocates its temporary source permutation.
    try std.testing.expectEqual(@as(usize, 3), identity_counter.alloc_index);
    try core.oom.checkAllocationFailures(
        std.testing.allocator,
        constructAndDiscardIdentity,
        .{source.len},
    );
}

test "input-order maps invert deterministic permutation families" {
    for (1..65) |count| {
        const source = try std.testing.allocator.alloc(core.ids.InputIndex, count);
        defer std.testing.allocator.free(source);
        const internal_values = try std.testing.allocator.alloc(u32, count);
        defer std.testing.allocator.free(internal_values);
        const output = try std.testing.allocator.alloc(u32, count);
        defer std.testing.allocator.free(output);

        const steps = [_]usize{ 1, count - 1 };
        const shifts = [_]usize{ 0, count / 2 };
        for (steps) |step| for (shifts) |shift| {
            for (source, internal_values, 0..) |*input_index, *value, internal_index| {
                input_index.* = @intCast((internal_index * step + shift) % count);
                value.* = @intCast(internal_index);
            }
            var order = try InputOrderMap.init(std.testing.allocator, source);
            defer order.deinit();
            try order.writeInputOrder(u32, internal_values, output);

            for (source, 0..) |input_index, internal_index| {
                try std.testing.expectEqual(@as(u32, @intCast(internal_index)), output[input_index]);
                try std.testing.expectEqual(@as(u32, @intCast(internal_index)), order.input_to_internal[input_index]);
            }
        };
    }
}
