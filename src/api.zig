const std = @import("std");
const core = @import("core.zig");

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
pub const Options = struct {
    precision: f32 = Precision.standard,
    score_residue_interactions: bool = true,
    treat_nonterminal_bonds_to_metal_as_zero_order: bool = true,
    even_angles: bool = false,
    skip_minimization: bool = false,
    force_open_macrocycles: bool = false,
    constrain_all_atoms: bool = false,
    build_from_fragments: bool = false,
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
        if (self.atoms.len > std.math.maxInt(u32) or
            self.bonds.len > std.math.maxInt(u32) or
            self.residues.len > std.math.maxInt(u32) or
            self.residue_interactions.len > std.math.maxInt(u32) or
            self.extra_bonds.len > std.math.maxInt(u32))
        {
            return error.TooManyItems;
        }
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
    try std.testing.expect(!options.build_from_fragments);
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
