//! The C ABI's data types, and nothing else.
//!
//! This file is the `c_abi` layer's root and defines no symbol that occupies
//! a linker name. The public entry points the DTOs describe live in
//! src/c_abi/exports.zig, in a module of their own, so that a binary can
//! import these declarations while linking somebody else's implementation of
//! coordgen_generate - which is precisely what the conformance corpus runner
//! does against the pinned C++ oracle. See that file's header for the full
//! reason and for cgz-r28.

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
    /// Formerly `build_from_fragments`; renamed by cgz-r25. The slot accepts
    /// exactly one value because upstream has no option here to gate (see
    /// `coordgen_options.reserved` in `include/coordgen_abi.h`). This struct
    /// is `extern` and bound to that frozen layout, so the slot stays at
    /// offset 28 under its honest name. Must be zero on input.
    reserved: u32 = 0,
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
    try expectLayout(core.math.Vec2, 8, 4);
    try expectOffset(core.math.Vec2, "x", 0);
    try expectOffset(core.math.Vec2, "y", 4);
    try expectLayout(core.math.Vec3, 12, 4);
    try expectOffset(core.math.Vec3, "x", 0);
    try expectOffset(core.math.Vec3, "y", 4);
    try expectOffset(core.math.Vec3, "z", 8);

    try expectLayout(StringView, 16, 8);
    try expectOffset(StringView, "ptr", 0);
    try expectOffset(StringView, "len", 8);
    try expectOffset(StringView, "reserved", 12);
    try expectLayout(IndexSpan, 16, 8);
    try expectOffset(IndexSpan, "ptr", 0);
    try expectOffset(IndexSpan, "len", 8);
    try expectOffset(IndexSpan, "reserved", 12);

    try expectLayout(AtomInput, 52, 4);
    try expectOffset(AtomInput, "atomic_number", 0);
    try expectOffset(AtomInput, "formal_charge", 4);
    try expectOffset(AtomInput, "flags", 8);
    try expectOffset(AtomInput, "stereo", 12);
    try expectOffset(AtomInput, "stereo_looking_from", 16);
    try expectOffset(AtomInput, "stereo_atom_a", 20);
    try expectOffset(AtomInput, "stereo_atom_b", 24);
    try expectOffset(AtomInput, "template_coordinates", 28);
    try expectOffset(AtomInput, "coordinates_3d", 36);
    try expectOffset(AtomInput, "reserved", 48);

    try expectLayout(BondInput, 40, 4);
    try expectOffset(BondInput, "start", 0);
    try expectOffset(BondInput, "end", 4);
    try expectOffset(BondInput, "order", 8);
    try expectOffset(BondInput, "flags", 12);
    try expectOffset(BondInput, "stereo", 16);
    try expectOffset(BondInput, "stereo_atom_a", 20);
    try expectOffset(BondInput, "stereo_atom_b", 24);
    try expectOffset(BondInput, "display", 28);
    try expectOffset(BondInput, "crossing_penalty_multiplier", 32);
    try expectOffset(BondInput, "reserved", 36);

    try expectLayout(ResidueInput, 32, 8);
    try expectOffset(ResidueInput, "atom", 0);
    try expectOffset(ResidueInput, "residue_number", 4);
    try expectOffset(ResidueInput, "closest_ligand_atom", 8);
    try expectOffset(ResidueInput, "reserved", 12);
    try expectOffset(ResidueInput, "chain", 16);

    try expectLayout(ResidueInteractionInput, 56, 8);
    try expectOffset(ResidueInteractionInput, "start", 0);
    try expectOffset(ResidueInteractionInput, "end", 4);
    try expectOffset(ResidueInteractionInput, "other_start_atoms", 8);
    try expectOffset(ResidueInteractionInput, "other_end_atoms", 24);
    try expectOffset(ResidueInteractionInput, "crossing_penalty_multiplier", 40);
    try expectOffset(ResidueInteractionInput, "reserved", 44);

    inline for (.{ AtomSpan, BondSpan, ResidueSpan, ResidueInteractionSpan, Vec2Span, U32Span }) |Span| {
        try expectLayout(Span, 16, 8);
        try expectOffset(Span, "ptr", 0);
        try expectOffset(Span, "len", 8);
        try expectOffset(Span, "reserved", 12);
    }

    try expectLayout(Options, 56, 8);
    try expectOffset(Options, "precision", 0);
    try expectOffset(Options, "score_residue_interactions", 4);
    try expectOffset(Options, "treat_nonterminal_bonds_to_metal_as_zero_order", 8);
    try expectOffset(Options, "even_angles", 12);
    try expectOffset(Options, "skip_minimization", 16);
    try expectOffset(Options, "force_open_macrocycles", 20);
    try expectOffset(Options, "constrain_all_atoms", 24);
    try expectOffset(Options, "reserved", 28);
    try expectOffset(Options, "debug_coordinates", 32);
    try expectOffset(Options, "load_templates", 36);
    try expectOffset(Options, "template_directory", 40);

    try expectLayout(Input, 136, 8);
    try expectOffset(Input, "options", 0);
    try expectOffset(Input, "atoms", 56);
    try expectOffset(Input, "bonds", 72);
    try expectOffset(Input, "residues", 88);
    try expectOffset(Input, "residue_interactions", 104);
    try expectOffset(Input, "extra_bonds", 120);

    try expectLayout(Result, 112, 8);
    try expectOffset(Result, "coordinates", 0);
    try expectOffset(Result, "input_to_internal", 16);
    try expectOffset(Result, "internal_to_input", 32);
    try expectOffset(Result, "effective_bond_orders", 48);
    try expectOffset(Result, "bond_displays", 64);
    try expectOffset(Result, "atom_stereo", 80);
    try expectOffset(Result, "clean_pose", 96);
    try expectOffset(Result, "reserved", 100);
    try expectOffset(Result, "owner", 104);
}

fn expectLayout(comptime T: type, size: usize, alignment: usize) !void {
    try std.testing.expectEqual(size, @sizeOf(T));
    try std.testing.expectEqual(alignment, @alignOf(T));
}

fn expectOffset(comptime T: type, comptime field: []const u8, offset: usize) !void {
    try std.testing.expectEqual(offset, @offsetOf(T, field));
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
    try std.testing.expectEqual(@as(u32, 0), options.reserved);
    try std.testing.expectEqual(@as(u32, 0), options.debug_coordinates);
    try std.testing.expectEqual(@as(u32, 1), options.load_templates);
    try std.testing.expect(options.template_directory.ptr == null);
}
