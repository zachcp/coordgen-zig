const std = @import("std");
const core = @import("core");

pub const abi_version: u32 = 1;
pub const bond_length: f32 = core.math.bond_length;

pub const AtomFlags = struct {
    pub const has_template_coordinates: u32 = 1 << 0;
    pub const has_coordinates_3d: u32 = 1 << 1;
    pub const fixed: u32 = 1 << 2;
    pub const constrained: u32 = 1 << 3;
    pub const hidden: u32 = 1 << 4;
    pub const known_mask: u32 = has_template_coordinates |
        has_coordinates_3d | fixed | constrained | hidden;
};

pub const BondFlags = struct {
    pub const skip: u32 = 1 << 0;
    pub const known_mask: u32 = skip;
};

pub const StringView = extern struct {
    ptr: ?[*]const u8 = null,
    len: u32 = 0,
    reserved: u32 = 0,
};

pub const IndexSpan = extern struct {
    ptr: ?[*]const u32 = null,
    len: u32 = 0,
    reserved: u32 = 0,
};

pub const AtomInput = extern struct {
    atomic_number: u32 = 6,
    formal_charge: i32 = 0,
    flags: u32 = 0,
    stereo: u32 = @backingInt(core.chemistry.AtomStereo.unspecified),
    stereo_looking_from: u32 = core.ids.invalid_input_index,
    stereo_atom_a: u32 = core.ids.invalid_input_index,
    stereo_atom_b: u32 = core.ids.invalid_input_index,
    template_coordinates: core.math.Vec2 = .{},
    coordinates_3d: core.math.Vec3 = .{},
    reserved: u32 = 0,
};

pub const BondInput = extern struct {
    start: u32,
    end: u32,
    order: u32 = @backingInt(core.chemistry.BondOrder.single),
    flags: u32 = 0,
    stereo: u32 = @backingInt(core.chemistry.BondStereo.unspecified),
    stereo_atom_a: u32 = core.ids.invalid_input_index,
    stereo_atom_b: u32 = core.ids.invalid_input_index,
    display: u32 = @backingInt(core.chemistry.BondDisplay.none),
    crossing_penalty_multiplier: f32 = 1,
    reserved: u32 = 0,
};

pub const ResidueInput = extern struct {
    atom: u32,
    residue_number: i32 = 0,
    closest_ligand_atom: u32 = core.ids.invalid_input_index,
    reserved: u32 = 0,
    chain: StringView = .{},
};

pub const ResidueInteractionInput = extern struct {
    start: u32,
    end: u32,
    other_start_atoms: IndexSpan = .{},
    other_end_atoms: IndexSpan = .{},
    crossing_penalty_multiplier: f32 = 1,
    reserved: [3]u32 = .{ 0, 0, 0 },
};

fn ConstSpan(comptime T: type) type {
    return extern struct {
        ptr: ?[*]const T = null,
        len: u32 = 0,
        reserved: u32 = 0,
    };
}

fn MutSpan(comptime T: type) type {
    return extern struct {
        ptr: ?[*]T = null,
        len: u32 = 0,
        reserved: u32 = 0,
    };
}

pub const AtomSpan = ConstSpan(AtomInput);
pub const BondSpan = ConstSpan(BondInput);
pub const ResidueSpan = ConstSpan(ResidueInput);
pub const ResidueInteractionSpan = ConstSpan(ResidueInteractionInput);

pub const Options = extern struct {
    precision: f32 = 1,
    score_residue_interactions: u32 = 1,
    treat_nonterminal_bonds_to_metal_as_zero_order: u32 = 1,
    even_angles: u32 = 0,
    skip_minimization: u32 = 0,
    force_open_macrocycles: u32 = 0,
    constrain_all_atoms: u32 = 0,
    build_from_fragments: u32 = 0,
    debug_coordinates: u32 = 0,
    load_templates: u32 = 1,
    template_directory: StringView = .{},
};

pub const Input = extern struct {
    options: Options = .{},
    atoms: AtomSpan = .{},
    bonds: BondSpan = .{},
    residues: ResidueSpan = .{},
    residue_interactions: ResidueInteractionSpan = .{},
    extra_bonds: BondSpan = .{},
};

pub const Vec2Span = MutSpan(core.math.Vec2);
pub const U32Span = MutSpan(u32);

/// On success, every span and `owner` are owned by the result and released by
/// coordgen_result_free(). On failure, the implementation leaves a zero value.
pub const Result = extern struct {
    coordinates: Vec2Span = .{},
    input_to_internal: U32Span = .{},
    internal_to_input: U32Span = .{},
    effective_bond_orders: U32Span = .{},
    bond_displays: U32Span = .{},
    atom_stereo: U32Span = .{},
    clean_pose: u32 = 0,
    reserved: u32 = 0,
    owner: ?*anyopaque = null,
};

test "C DTO widths, alignments, and offsets are frozen on supported 64-bit targets" {
    try std.testing.expectEqual(@as(usize, 16), @sizeOf(StringView));
    try std.testing.expectEqual(@as(usize, 52), @sizeOf(AtomInput));
    try std.testing.expectEqual(@as(usize, 40), @sizeOf(BondInput));
    try std.testing.expectEqual(@as(usize, 32), @sizeOf(ResidueInput));
    try std.testing.expectEqual(@as(usize, 56), @sizeOf(ResidueInteractionInput));
    try std.testing.expectEqual(@as(usize, 56), @sizeOf(Options));
    try std.testing.expectEqual(@as(usize, 136), @sizeOf(Input));
    try std.testing.expectEqual(@as(usize, 112), @sizeOf(Result));
    try std.testing.expectEqual(@as(usize, 8), @alignOf(Input));
    try std.testing.expectEqual(@as(usize, 28), @offsetOf(AtomInput, "template_coordinates"));
    try std.testing.expectEqual(@as(usize, 48), @offsetOf(AtomInput, "reserved"));
    try std.testing.expectEqual(@as(usize, 40), @offsetOf(Options, "template_directory"));
    try std.testing.expectEqual(@as(usize, 104), @offsetOf(Result, "owner"));
}

test "ABI numeric values match conserved chemistry types" {
    try std.testing.expectEqual(@as(u32, 0), @backingInt(core.chemistry.BondOrder.zero));
    try std.testing.expectEqual(@as(u32, 4), @backingInt(core.chemistry.BondDisplay.hashed_reverse));
    try std.testing.expectEqual(@as(u32, 12), @backingInt(core.errors.ErrorCode.internal));
}

test "C option DTO defaults match the safe Zig option contract" {
    const options = Options{};
    try std.testing.expectEqual(@as(f32, 1), options.precision);
    try std.testing.expectEqual(@as(u32, 1), options.score_residue_interactions);
    try std.testing.expectEqual(@as(u32, 1), options.treat_nonterminal_bonds_to_metal_as_zero_order);
    try std.testing.expectEqual(@as(u32, 0), options.even_angles);
    try std.testing.expectEqual(@as(u32, 0), options.skip_minimization);
    try std.testing.expectEqual(@as(u32, 0), options.force_open_macrocycles);
    try std.testing.expectEqual(@as(u32, 0), options.constrain_all_atoms);
    try std.testing.expectEqual(@as(u32, 0), options.build_from_fragments);
    try std.testing.expectEqual(@as(u32, 0), options.debug_coordinates);
    try std.testing.expectEqual(@as(u32, 1), options.load_templates);
    try std.testing.expect(options.template_directory.ptr == null);
}
