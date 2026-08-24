const std = @import("std");
const core = @import("core");
const generator = @import("generator");

pub const InputIndex = core.ids.InputIndex;
pub const invalid_input_index = core.ids.invalid_input_index;
pub const AtomicNumber = core.chemistry.AtomicNumber;
pub const FormalCharge = core.chemistry.FormalCharge;
pub const BondOrder = core.chemistry.BondOrder;
pub const AtomStereo = core.chemistry.AtomStereo;
pub const BondStereo = core.chemistry.BondStereo;
pub const BondDisplay = core.chemistry.BondDisplay;
pub const Vec2 = core.math.Vec2;
pub const Vec3 = core.math.Vec3;
pub const bond_length = core.math.bond_length;

pub const Precision = struct {
    pub const quick: f32 = 0.2;
    pub const standard: f32 = 1.0;
    pub const best: f32 = 3.0;
};

pub const AtomStereoInput = struct {
    value: AtomStereo = .unspecified,
    looking_from: ?InputIndex = null,
    atom_a: ?InputIndex = null,
    atom_b: ?InputIndex = null,
};

pub const AtomInput = struct {
    atomic_number: AtomicNumber = .carbon,
    formal_charge: FormalCharge = 0,
    fixed: bool = false,
    constrained: bool = false,
    hidden: bool = false,
    template_coordinates: ?Vec2 = null,
    coordinates_3d: ?Vec3 = null,
    stereo: AtomStereoInput = .{},
};

pub const BondStereoInput = struct {
    value: BondStereo = .unspecified,
    atom_a: ?InputIndex = null,
    atom_b: ?InputIndex = null,
};

pub const BondInput = struct {
    start: InputIndex,
    end: InputIndex,
    order: BondOrder = .single,
    skip: bool = false,
    stereo: BondStereoInput = .{},
    display: BondDisplay = .none,
    crossing_penalty_multiplier: f32 = 1.0,
};

pub const ResidueInput = struct {
    atom: InputIndex,
    chain: []const u8 = "",
    residue_number: i32 = 0,
    closest_ligand_atom: ?InputIndex = null,
};

pub const ResidueInteractionInput = struct {
    start: InputIndex,
    end: InputIndex,
    other_start_atoms: []const InputIndex = &.{},
    other_end_atoms: []const InputIndex = &.{},
    crossing_penalty_multiplier: f32 = 1.0,
};

/// All upstream generation controls are represented in a single-shot request.
/// Per-atom constrain/fix vectors are normalized into AtomInput flags;
/// addExtraBond is normalized into Input.extra_bonds.
///
/// There is no `build_from_fragments` field. cgz-r25 removed it rather than
/// mirror `c_abi.Options.reserved` here: upstream's buildFromFragments(bool)
/// is an imperative pipeline step runGenerateCoordinates() already performs
/// unconditionally (cgz-7v2.8), never a stored option, so the field could
/// accept only its default. `c_abi.Options` keeps the slot because it is
/// `extern` and bound by coordgen_abi.h's frozen byte layout; this struct is
/// plain Zig with no layout obligation, so the slot has nothing holding it in
/// place and is simply gone.
pub const Options = struct {
    precision: f32 = Precision.standard,
    score_residue_interactions: bool = true,
    treat_nonterminal_bonds_to_metal_as_zero_order: bool = true,
    even_angles: bool = false,
    skip_minimization: bool = false,
    force_open_macrocycles: bool = false,
    constrain_all_atoms: bool = false,
    debug_coordinates: bool = false,
    load_templates: bool = true,
    /// Borrowed for the duration of generate(); null selects built-in data.
    template_directory: ?[]const u8 = null,
};

/// Input arrays and strings are immutable and borrowed only for generate().
/// A one-shot graph makes the upstream bond/finalize/stereo call-order hazard
/// unrepresentable in both Zig and C.
pub const Input = struct {
    atoms: []const AtomInput,
    bonds: []const BondInput,
    residues: []const ResidueInput = &.{},
    residue_interactions: []const ResidueInteractionInput = &.{},
    extra_bonds: []const BondInput = &.{},
    options: Options = .{},

    pub fn validate(self: Input) core.errors.Error!void {
        if (self.atoms.len == 0) return error.EmptyGraph;
        try validateItemCounts(
            self.atoms.len,
            self.bonds.len,
            self.residues.len,
            self.residue_interactions.len,
            self.extra_bonds.len,
        );
        if (!std.math.isFinite(self.options.precision) or self.options.precision <= 0) {
            return error.InvalidOption;
        }

        for (self.atoms) |atom| {
            const number = @backingInt(atom.atomic_number);
            if (number < 1 or number > 118) return error.InvalidAtomicNumber;
            if (atom.template_coordinates) |position| {
                if (!position.isFinite()) return error.InvalidCoordinate;
            }
            if (atom.coordinates_3d) |position| {
                if (!position.isFinite()) return error.InvalidCoordinate;
            }
            try validateAtomStereo(atom.stereo, self.atoms.len);
        }
        for (self.bonds) |bond| try validateBond(bond, self.atoms.len);
        for (self.extra_bonds) |bond| try validateBond(bond, self.atoms.len);
        for (self.residues) |residue| {
            try validateIndex(residue.atom, self.atoms.len);
            if (residue.closest_ligand_atom) |index| try validateIndex(index, self.atoms.len);
        }
        for (self.residue_interactions) |interaction| {
            try validateIndex(interaction.start, self.atoms.len);
            try validateIndex(interaction.end, self.atoms.len);
            if (!std.math.isFinite(interaction.crossing_penalty_multiplier) or
                interaction.crossing_penalty_multiplier < 0)
            {
                return error.InvalidOption;
            }
            for (interaction.other_start_atoms) |index| try validateIndex(index, self.atoms.len);
            for (interaction.other_end_atoms) |index| try validateIndex(index, self.atoms.len);
        }
    }
};

fn validateItemCounts(
    atom_count: usize,
    bond_count: usize,
    residue_count: usize,
    interaction_count: usize,
    extra_bond_count: usize,
) core.errors.Error!void {
    if (atom_count > std.math.maxInt(u32) or
        bond_count > std.math.maxInt(u32) or
        residue_count > std.math.maxInt(u32) or
        interaction_count > std.math.maxInt(u32) or
        extra_bond_count > std.math.maxInt(u32))
    {
        return error.TooManyItems;
    }
}

fn validateIndex(index: InputIndex, count: usize) core.errors.Error!void {
    if (index == invalid_input_index or index >= count) return error.InvalidAtomIndex;
}

fn validateAtomStereo(stereo: AtomStereoInput, atom_count: usize) core.errors.Error!void {
    switch (stereo.value) {
        .unspecified => if (stereo.looking_from != null or stereo.atom_a != null or stereo.atom_b != null) {
            return error.InvalidStereo;
        },
        .clockwise, .counter_clockwise => {
            try validateIndex(stereo.looking_from orelse return error.InvalidStereo, atom_count);
            try validateIndex(stereo.atom_a orelse return error.InvalidStereo, atom_count);
            try validateIndex(stereo.atom_b orelse return error.InvalidStereo, atom_count);
        },
        .r, .s => if (stereo.looking_from != null or stereo.atom_a != null or stereo.atom_b != null) {
            return error.InvalidStereo;
        },
    }
}

fn validateBond(bond: BondInput, atom_count: usize) core.errors.Error!void {
    try validateIndex(bond.start, atom_count);
    try validateIndex(bond.end, atom_count);
    if (bond.start == bond.end) return error.InvalidAtomIndex;
    if (core.chemistry.BondOrder.fromInt(@backingInt(bond.order)) == null) {
        return error.InvalidBondOrder;
    }
    if (!std.math.isFinite(bond.crossing_penalty_multiplier) or
        bond.crossing_penalty_multiplier < 0)
    {
        return error.InvalidOption;
    }
    switch (bond.stereo.value) {
        .unspecified => if (bond.stereo.atom_a != null or bond.stereo.atom_b != null) {
            return error.InvalidStereo;
        },
        .cis, .trans => {
            try validateIndex(bond.stereo.atom_a orelse return error.InvalidStereo, atom_count);
            try validateIndex(bond.stereo.atom_b orelse return error.InvalidStereo, atom_count);
        },
        .z, .e => if (bond.stereo.atom_a != null or bond.stereo.atom_b != null) {
            return error.InvalidStereo;
        },
    }
}

/// Allocator-owned result. Every slice is indexed by caller input order except
/// the two explicitly named permutation maps. `deinit` must be called exactly
/// once with the allocator retained in this value.
pub const Result = struct {
    allocator: std.mem.Allocator,
    coordinates: []Vec2,
    input_to_internal: []u32,
    internal_to_input: []InputIndex,
    effective_bond_orders: []BondOrder,
    bond_displays: []BondDisplay,
    atom_stereo: []AtomStereo,
    clean_pose: bool = false,

    pub fn init(allocator: std.mem.Allocator, atom_count: usize, bond_count: usize) !Result {
        const coordinates = try allocator.alloc(Vec2, atom_count);
        errdefer allocator.free(coordinates);
        const input_to_internal = try allocator.alloc(u32, atom_count);
        errdefer allocator.free(input_to_internal);
        const internal_to_input = try allocator.alloc(InputIndex, atom_count);
        errdefer allocator.free(internal_to_input);
        const effective_bond_orders = try allocator.alloc(BondOrder, bond_count);
        errdefer allocator.free(effective_bond_orders);
        const bond_displays = try allocator.alloc(BondDisplay, bond_count);
        errdefer allocator.free(bond_displays);
        const atom_stereo = try allocator.alloc(AtomStereo, atom_count);
        errdefer allocator.free(atom_stereo);
        return .{
            .allocator = allocator,
            .coordinates = coordinates,
            .input_to_internal = input_to_internal,
            .internal_to_input = internal_to_input,
            .effective_bond_orders = effective_bond_orders,
            .bond_displays = bond_displays,
            .atom_stereo = atom_stereo,
        };
    }

    pub fn deinit(self: *Result) void {
        self.allocator.free(self.atom_stereo);
        self.allocator.free(self.bond_displays);
        self.allocator.free(self.effective_bond_orders);
        self.allocator.free(self.internal_to_input);
        self.allocator.free(self.input_to_internal);
        self.allocator.free(self.coordinates);
        self.* = undefined;
    }
};

/// Generate 2D coordinates for `input`. The caller owns the returned value and
/// must call `deinit` on it exactly once with the same allocator.
///
/// `input` is borrowed for the duration of the call and is never mutated:
/// upstream's `initialize()` takes ownership of caller atoms and rewrites
/// nonterminal metal bond orders in place, and this port surfaces that
/// normalization as `Result.effective_bond_orders` instead (cgz-r11).
///
/// Coverage is deliberately partial and deliberately loud. An input whose
/// domain no native phase owns yet returns `error.Unsupported` rather than a
/// successful empty result; the exact rejected set is
/// `generator.rejectOutOfScope`, and every rejection is a domain, never a
/// size or a difficulty. Supported today: ordinary bonds, zero-order and
/// residue-interaction proximity, rings and macrocycles, templates, residues,
/// and the discrete pose search. Rejected today: extra bonds, hidden, fixed,
/// constrained, template-positioned, or 3D-seeded atoms, caller-specified
/// atom and bond stereo, skipped bonds, caller bond displays, and the
/// `even_angles`, `constrain_all_atoms`, `debug_coordinates`, and
/// `template_directory` options. The C ABI's `reserved` slot is rejected at
/// that boundary instead, since cgz-r25 removed it from this option table.
///
/// On any error the result is released before returning, so a failed call
/// leaks nothing and leaves the caller nothing to free.
pub fn generate(allocator: std.mem.Allocator, input: Input) core.errors.Error!Result {
    try input.validate();
    var result = Result.init(allocator, input.atoms.len, input.bonds.len) catch {
        return error.OutOfMemory;
    };
    errdefer result.deinit();
    result.clean_pose = try generator.generateInto(allocator, input, .{
        .coordinates = result.coordinates,
        .input_to_internal = result.input_to_internal,
        .internal_to_input = result.internal_to_input,
        .effective_bond_orders = result.effective_bond_orders,
        .bond_displays = result.bond_displays,
        .atom_stereo = result.atom_stereo,
    });
    return result;
}

test "public generation emits caller-ordered coordinates and every normalization output" {
    const atoms = [_]AtomInput{
        .{ .atomic_number = .carbon },
        .{ .atomic_number = .oxygen },
        .{ .atomic_number = .nitrogen },
    };
    const bonds = [_]BondInput{
        .{ .start = 0, .end = 1, .order = .single },
        .{ .start = 1, .end = 2, .order = .double },
    };
    const atoms_before = atoms;
    const bonds_before = bonds;
    const input = Input{ .atoms = &atoms, .bonds = &bonds };

    var result = try generate(std.testing.allocator, input);
    defer result.deinit();

    // Borrowed, never mutated: the upstream hazard cgz-r11 records.
    try std.testing.expectEqualDeep(atoms_before, atoms);
    try std.testing.expectEqualDeep(bonds_before, bonds);

    try std.testing.expectEqual(atoms.len, result.coordinates.len);
    try std.testing.expectEqual(bonds.len, result.effective_bond_orders.len);
    try std.testing.expectEqual(bonds.len, result.bond_displays.len);
    try std.testing.expectEqual(atoms.len, result.atom_stereo.len);

    for (result.coordinates) |position| try std.testing.expect(position.isFinite());
    for (bonds) |bond| {
        const delta_x = result.coordinates[bond.start].x - result.coordinates[bond.end].x;
        const delta_y = result.coordinates[bond.start].y - result.coordinates[bond.end].y;
        try std.testing.expectApproxEqAbs(bond_length, @sqrt(delta_x * delta_x + delta_y * delta_y), 0.001);
    }
    // A successful call that returned zeroed storage would pass every length
    // assertion above; these are the facts that separate the two.
    for (bonds, result.effective_bond_orders) |bond, effective| {
        try std.testing.expectEqual(bond.order, effective);
    }
    for (result.bond_displays) |display| try std.testing.expectEqual(BondDisplay.none, display);
    for (result.atom_stereo) |stereo| try std.testing.expectEqual(AtomStereo.unspecified, stereo);
    for (result.internal_to_input, 0..) |input_index, internal_index| {
        try std.testing.expectEqual(@as(u32, @intCast(internal_index)), result.input_to_internal[input_index]);
    }
}

test "public generation surfaces the metal zero-order rewrite instead of mutating the caller" {
    const atoms = [_]AtomInput{
        .{ .atomic_number = .carbon },
        .{ .atomic_number = .carbon },
        .{ .atomic_number = .iron },
        .{ .atomic_number = .carbon },
        .{ .atomic_number = .carbon },
    };
    // Only the two bonds at the iron have a nonterminal atom at both ends,
    // which is the condition the rewrite applies to.
    const bonds = [_]BondInput{
        .{ .start = 0, .end = 1, .order = .single },
        .{ .start = 1, .end = 2, .order = .single },
        .{ .start = 2, .end = 3, .order = .single },
        .{ .start = 3, .end = 4, .order = .single },
    };
    const expected = [_]BondOrder{ .single, .zero, .zero, .single };
    const input = Input{ .atoms = &atoms, .bonds = &bonds };

    var rewritten = try generate(std.testing.allocator, input);
    defer rewritten.deinit();
    try std.testing.expectEqualSlices(BondOrder, &expected, rewritten.effective_bond_orders);
    for (bonds) |bond| try std.testing.expectEqual(BondOrder.single, bond.order);

    var preserved = try generate(std.testing.allocator, .{
        .atoms = &atoms,
        .bonds = &bonds,
        .options = .{ .treat_nonterminal_bonds_to_metal_as_zero_order = false },
    });
    defer preserved.deinit();
    for (preserved.effective_bond_orders) |effective| {
        try std.testing.expectEqual(BondOrder.single, effective);
    }
}

test "public generation rejects an unowned domain instead of reporting empty success" {
    const atoms = [_]AtomInput{ .{}, .{}, .{} };
    const bonds = [_]BondInput{ .{ .start = 0, .end = 1 }, .{ .start = 1, .end = 2 } };
    const extra = [_]BondInput{.{ .start = 0, .end = 2 }};
    try std.testing.expectError(error.Unsupported, generate(std.testing.allocator, .{
        .atoms = &atoms,
        .bonds = &bonds,
        .extra_bonds = &extra,
    }));

    const hidden = [_]AtomInput{ .{}, .{}, .{ .hidden = true } };
    try std.testing.expectError(error.Unsupported, generate(std.testing.allocator, .{
        .atoms = &hidden,
        .bonds = &bonds,
    }));

    // Validation errors still precede domain rejection.
    try std.testing.expectError(error.EmptyGraph, generate(std.testing.allocator, .{
        .atoms = &.{},
        .bonds = &.{},
    }));
}

fn generateAndDiscard(allocator: std.mem.Allocator, input: Input) !void {
    var result = try generate(allocator, input);
    result.deinit();
}

test "public generation reports and cleans up every allocation failure" {
    const atoms = [_]AtomInput{ .{}, .{}, .{} };
    const bonds = [_]BondInput{ .{ .start = 0, .end = 1 }, .{ .start = 1, .end = 2 } };
    const input = Input{ .atoms = &atoms, .bonds = &bonds };
    // checkAllAllocationFailures requires a deterministic allocation count and
    // reports a first-call warm-up as nondeterminism. Discharge it here, and
    // prove it is a warm-up rather than a leak by requiring the two following
    // measured runs to agree exactly.
    try generateAndDiscard(std.testing.allocator, input);
    var baseline = std.testing.FailingAllocator.init(std.testing.allocator, .{});
    try generateAndDiscard(baseline.allocator(), input);
    var repeated = std.testing.FailingAllocator.init(std.testing.allocator, .{});
    try generateAndDiscard(repeated.allocator(), input);
    try std.testing.expectEqual(baseline.alloc_index, repeated.alloc_index);
    try std.testing.expectEqual(baseline.allocated_bytes, baseline.freed_bytes);

    try core.oom.checkAllocationFailures(
        std.testing.allocator,
        generateAndDiscard,
        .{input},
    );
}

test "bond length and named precision values are conserved" {
    try std.testing.expectEqual(@as(f32, 50), bond_length);
    try std.testing.expectEqual(@as(f32, 0.2), Precision.quick);
    try std.testing.expectEqual(@as(f32, 1.0), Precision.standard);
    try std.testing.expectEqual(@as(f32, 3.0), Precision.best);
}

test "complete option table has pinned defaults" {
    const options = Options{};
    try std.testing.expectEqual(Precision.standard, options.precision);
    try std.testing.expect(options.score_residue_interactions);
    try std.testing.expect(options.treat_nonterminal_bonds_to_metal_as_zero_order);
    try std.testing.expect(!options.even_angles);
    try std.testing.expect(!options.skip_minimization);
    try std.testing.expect(!options.force_open_macrocycles);
    try std.testing.expect(!options.constrain_all_atoms);
    try std.testing.expect(!options.debug_coordinates);
    try std.testing.expect(options.load_templates);
    try std.testing.expect(options.template_directory == null);
}

test "zero-order, hidden, fixed, constrained, templates, charge, and stereo validate" {
    const atoms = [_]AtomInput{
        .{
            .atomic_number = .carbon,
            .formal_charge = -1,
            .fixed = true,
            .constrained = true,
            .hidden = true,
            .template_coordinates = .{ .x = 1, .y = 2 },
        },
        .{ .atomic_number = .oxygen },
        .{ .atomic_number = .hydrogen },
        .{ .atomic_number = .hydrogen },
    };
    const bonds = [_]BondInput{.{
        .start = 0,
        .end = 1,
        .order = .zero,
        .skip = true,
        .stereo = .{ .value = .cis, .atom_a = 2, .atom_b = 3 },
        .display = .hashed_reverse,
    }};
    try (Input{ .atoms = &atoms, .bonds = &bonds }).validate();
}

test "result allocation has explicit leak-free ownership" {
    var result = try Result.init(std.testing.allocator, 3, 2);
    result.deinit();
}

fn allocateAndDiscardResult(
    allocator: std.mem.Allocator,
    atom_count: usize,
    bond_count: usize,
) !void {
    var result = try Result.init(allocator, atom_count, bond_count);
    result.deinit();
}

test "result reports and cleans up every allocation failure" {
    var counting_allocator = std.testing.FailingAllocator.init(std.testing.allocator, .{});
    try allocateAndDiscardResult(counting_allocator.allocator(), 3, 2);
    // One independently-fallible allocation for each owned output slice. A
    // decrease means an output stopped participating in the OOM sweep.
    try std.testing.expectEqual(@as(usize, 6), counting_allocator.alloc_index);
    try core.oom.checkAllocationFailures(
        std.testing.allocator,
        allocateAndDiscardResult,
        .{ 3, 2 },
    );
}

test "caller-controlled invalid states return errors" {
    try std.testing.expectError(error.EmptyGraph, (Input{
        .atoms = &.{},
        .bonds = &.{},
    }).validate());

    const atoms = [_]AtomInput{
        .{ .atomic_number = .carbon },
        .{ .atomic_number = .oxygen },
    };
    try std.testing.expectError(error.InvalidAtomIndex, (Input{
        .atoms = &atoms,
        .bonds = &.{.{ .start = 0, .end = 2 }},
    }).validate());
    try std.testing.expectError(error.InvalidOption, (Input{
        .atoms = &atoms,
        .bonds = &.{.{ .start = 0, .end = 1 }},
        .options = .{ .precision = std.math.nan(f32) },
    }).validate());
}

test "all safe API malformed fields return their documented errors" {
    const atoms = [_]AtomInput{
        .{ .atomic_number = .carbon },
        .{ .atomic_number = .oxygen },
        .{ .atomic_number = .hydrogen },
        .{ .atomic_number = .hydrogen },
    };

    try std.testing.expectError(error.InvalidAtomicNumber, (Input{
        .atoms = &.{.{ .atomic_number = .virtual }},
        .bonds = &.{},
    }).validate());
    try std.testing.expectError(error.InvalidCoordinate, (Input{
        .atoms = &.{.{ .template_coordinates = .{ .x = std.math.nan(f32), .y = 0 } }},
        .bonds = &.{},
    }).validate());
    try std.testing.expectError(error.InvalidCoordinate, (Input{
        .atoms = &.{.{ .coordinates_3d = .{ .x = 0, .y = 0, .z = std.math.inf(f32) } }},
        .bonds = &.{},
    }).validate());

    const invalid_atom_stereo = [_]struct { stereo: AtomStereoInput, expected: core.errors.Error }{
        .{ .stereo = .{ .looking_from = 0 }, .expected = error.InvalidStereo },
        .{ .stereo = .{ .value = .clockwise }, .expected = error.InvalidStereo },
        .{
            .stereo = .{ .value = .clockwise, .looking_from = 0, .atom_a = 1, .atom_b = 9 },
            .expected = error.InvalidAtomIndex,
        },
        .{ .stereo = .{ .value = .r, .atom_a = 1 }, .expected = error.InvalidStereo },
    };
    for (invalid_atom_stereo) |case| {
        var changed = atoms;
        changed[0].stereo = case.stereo;
        try std.testing.expectError(case.expected, (Input{
            .atoms = &changed,
            .bonds = &.{},
        }).validate());
    }

    const invalid_bonds = [_]struct { bond: BondInput, expected: core.errors.Error }{
        .{ .bond = .{ .start = 0, .end = 0 }, .expected = error.InvalidAtomIndex },
        .{ .bond = .{ .start = 0, .end = 9 }, .expected = error.InvalidAtomIndex },
        .{
            .bond = .{ .start = 0, .end = 1, .crossing_penalty_multiplier = -1 },
            .expected = error.InvalidOption,
        },
        .{
            .bond = .{ .start = 0, .end = 1, .crossing_penalty_multiplier = std.math.inf(f32) },
            .expected = error.InvalidOption,
        },
        .{
            .bond = .{ .start = 0, .end = 1, .stereo = .{ .atom_a = 2 } },
            .expected = error.InvalidStereo,
        },
        .{
            .bond = .{ .start = 0, .end = 1, .stereo = .{ .value = .cis, .atom_a = 2 } },
            .expected = error.InvalidStereo,
        },
        .{
            .bond = .{ .start = 0, .end = 1, .stereo = .{ .value = .trans, .atom_a = 2, .atom_b = 9 } },
            .expected = error.InvalidAtomIndex,
        },
        .{
            .bond = .{ .start = 0, .end = 1, .stereo = .{ .value = .z, .atom_a = 2 } },
            .expected = error.InvalidStereo,
        },
    };
    for (invalid_bonds) |case| {
        try std.testing.expectError(case.expected, (Input{
            .atoms = &atoms,
            .bonds = @as(*const [1]BondInput, &case.bond),
        }).validate());
    }

    try std.testing.expectError(error.InvalidAtomIndex, (Input{
        .atoms = &atoms,
        .bonds = &.{},
        .extra_bonds = &.{.{ .start = 0, .end = 9 }},
    }).validate());
    try std.testing.expectError(error.InvalidAtomIndex, (Input{
        .atoms = &atoms,
        .bonds = &.{},
        .residues = &.{.{ .atom = 9 }},
    }).validate());
    try std.testing.expectError(error.InvalidAtomIndex, (Input{
        .atoms = &atoms,
        .bonds = &.{},
        .residues = &.{.{ .atom = 0, .closest_ligand_atom = 9 }},
    }).validate());

    const invalid_interactions = [_]ResidueInteractionInput{
        .{ .start = 9, .end = 1 },
        .{ .start = 0, .end = 9 },
        .{ .start = 0, .end = 1, .crossing_penalty_multiplier = -1 },
        .{ .start = 0, .end = 1, .crossing_penalty_multiplier = std.math.nan(f32) },
        .{ .start = 0, .end = 1, .other_start_atoms = &.{9} },
        .{ .start = 0, .end = 1, .other_end_atoms = &.{9} },
    };
    for (invalid_interactions) |interaction| {
        const expected: core.errors.Error = if (!std.math.isFinite(interaction.crossing_penalty_multiplier) or
            interaction.crossing_penalty_multiplier < 0)
            error.InvalidOption
        else
            error.InvalidAtomIndex;
        try std.testing.expectError(expected, (Input{
            .atoms = &atoms,
            .bonds = &.{},
            .residue_interactions = @as(*const [1]ResidueInteractionInput, &interaction),
        }).validate());
    }

    for ([_]f32{ 0, -1, std.math.inf(f32), std.math.nan(f32) }) |precision| {
        try std.testing.expectError(error.InvalidOption, (Input{
            .atoms = &atoms,
            .bonds = &.{},
            .options = .{ .precision = precision },
        }).validate());
    }
}

test "every item-count field rejects values wider than the public index type" {
    if (@sizeOf(usize) <= @sizeOf(u32)) return;
    const too_many = @as(usize, std.math.maxInt(u32)) + 1;
    inline for (0..5) |position| {
        var counts = [_]usize{ 1, 0, 0, 0, 0 };
        counts[position] = too_many;
        try std.testing.expectError(
            error.TooManyItems,
            validateItemCounts(counts[0], counts[1], counts[2], counts[3], counts[4]),
        );
    }
}
