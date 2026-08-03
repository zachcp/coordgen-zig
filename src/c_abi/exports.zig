//! Sole definition site of the public C ABI entry points declared in
//! include/coordgen_abi.h.
//!
//! These live apart from the DTOs in src/c_abi_types.zig for a link-level
//! reason, not a stylistic one. `conformance/oracle_adapter.cpp` implements
//! the *same* symbol names on purpose, so that a C or C++ consumer of the
//! public header can be linked against the pinned upstream oracle without
//! being rewritten. Any binary that wants the oracle's implementation and
//! also wants the ABI's Zig DTO declarations - `tests/oracle_corpus_run.zig`
//! is exactly that - would link two definitions of coordgen_generate if the
//! DTO module carried the exports too. That is a hard `ld.lld: duplicate
//! symbol` error on Linux and silently tolerated by macOS's linker; see
//! cgz-r28.
//!
//! So the ownership rule is: the native library owns these names, and it is
//! the only Zig-side provider of them. Importing `c_abi` gets you types and
//! nothing that occupies a public symbol. Only src/coordgen.zig, the
//! installed library's root, pulls this module in - and it must do so with a
//! `comptime` reference, because a plain `pub const` @import does not force
//! export discovery in a non-test build (that is the cgz-r28 defect itself).

const std = @import("std");
const core = @import("core");
const c_abi = @import("c_abi");

const AtomFlags = c_abi.AtomFlags;
const BondFlags = c_abi.BondFlags;
const AtomInput = c_abi.AtomInput;
const BondInput = c_abi.BondInput;
const ResidueInput = c_abi.ResidueInput;
const ResidueInteractionInput = c_abi.ResidueInteractionInput;
const Input = c_abi.Input;
const Result = c_abi.Result;

fn asSlice(comptime T: type, span: anytype) []const T {
    const len: usize = @intCast(span.len);
    return if (span.ptr) |p| p[0..len] else &[_]T{};
}

fn isIndexPresent(value: u32) bool {
    return value != core.ids.invalid_input_index;
}

fn validateAtomStereoField(atom: AtomInput, atom_count: usize) ?core.errors.ErrorCode {
    const stereo = core.chemistry.AtomStereo.fromPublic(atom.stereo) orelse return .invalid_stereo;
    switch (stereo) {
        .unspecified, .r, .s => {
            if (isIndexPresent(atom.stereo_looking_from) or
                isIndexPresent(atom.stereo_atom_a) or
                isIndexPresent(atom.stereo_atom_b))
            {
                return .invalid_stereo;
            }
        },
        .clockwise, .counter_clockwise => {
            if (!isIndexPresent(atom.stereo_looking_from) or
                !isIndexPresent(atom.stereo_atom_a) or
                !isIndexPresent(atom.stereo_atom_b))
            {
                return .invalid_stereo;
            }
            if (atom.stereo_looking_from >= atom_count or
                atom.stereo_atom_a >= atom_count or
                atom.stereo_atom_b >= atom_count)
            {
                return .invalid_atom_index;
            }
        },
    }
    return null;
}

/// Every check here mirrors api.Input.validate()'s atom rules exactly, but
/// operates on the C sentinel/flag representation directly: this is the
/// c_abi layer's own job (FOUNDATION_CONTRACTS.md lists "unknown C flag bits
/// rejected" under the C ABI column specifically, since a safe Zig caller
/// cannot construct an unknown flag bit in the first place).
fn validateAtomField(atom: AtomInput, atom_count: usize) ?core.errors.ErrorCode {
    if (atom.flags & ~AtomFlags.known_mask != 0) return .invalid_option;
    if (core.chemistry.AtomicNumber.fromPublic(atom.atomic_number) == null) return .invalid_atomic_number;
    if (atom.flags & AtomFlags.has_template_coordinates != 0 and !atom.template_coordinates.isFinite()) {
        return .invalid_coordinate;
    }
    if (atom.flags & AtomFlags.has_coordinates_3d != 0 and !atom.coordinates_3d.isFinite()) {
        return .invalid_coordinate;
    }
    if (validateAtomStereoField(atom, atom_count)) |err| return err;
    return null;
}

fn validateBondStereoField(bond: BondInput, atom_count: usize) ?core.errors.ErrorCode {
    const stereo = core.chemistry.BondStereo.fromPublic(bond.stereo) orelse return .invalid_stereo;
    switch (stereo) {
        .unspecified, .z, .e => {
            if (isIndexPresent(bond.stereo_atom_a) or isIndexPresent(bond.stereo_atom_b)) {
                return .invalid_stereo;
            }
        },
        .cis, .trans => {
            if (!isIndexPresent(bond.stereo_atom_a) or !isIndexPresent(bond.stereo_atom_b)) {
                return .invalid_stereo;
            }
            if (bond.stereo_atom_a >= atom_count or bond.stereo_atom_b >= atom_count) {
                return .invalid_atom_index;
            }
        },
    }
    return null;
}

/// Mirrors api.zig's validateBond(); see validateAtomField for why this
/// lives here rather than reusing that function directly (it operates on
/// the safe Zig Input, not the borrowed C span).
fn validateBondField(bond: BondInput, atom_count: usize) ?core.errors.ErrorCode {
    if (bond.flags & ~BondFlags.known_mask != 0) return .invalid_option;
    if (!isIndexPresent(bond.start) or bond.start >= atom_count) return .invalid_atom_index;
    if (!isIndexPresent(bond.end) or bond.end >= atom_count) return .invalid_atom_index;
    if (bond.start == bond.end) return .invalid_atom_index;
    if (core.chemistry.BondOrder.fromInt(bond.order) == null) return .invalid_bond_order;
    if (!std.math.isFinite(bond.crossing_penalty_multiplier) or bond.crossing_penalty_multiplier < 0) {
        return .invalid_option;
    }
    if (validateBondStereoField(bond, atom_count)) |err| return err;
    return null;
}

fn validateResidueField(residue: ResidueInput, atom_count: usize) ?core.errors.ErrorCode {
    if (!isIndexPresent(residue.atom) or residue.atom >= atom_count) return .invalid_atom_index;
    if (isIndexPresent(residue.closest_ligand_atom) and residue.closest_ligand_atom >= atom_count) {
        return .invalid_atom_index;
    }
    return null;
}

fn validateResidueInteractionField(
    interaction: ResidueInteractionInput,
    atom_count: usize,
) ?core.errors.ErrorCode {
    if (!isIndexPresent(interaction.start) or interaction.start >= atom_count) return .invalid_atom_index;
    if (!isIndexPresent(interaction.end) or interaction.end >= atom_count) return .invalid_atom_index;
    if (!std.math.isFinite(interaction.crossing_penalty_multiplier) or
        interaction.crossing_penalty_multiplier < 0)
    {
        return .invalid_option;
    }
    for (asSlice(u32, interaction.other_start_atoms)) |index| {
        if (!isIndexPresent(index) or index >= atom_count) return .invalid_atom_index;
    }
    for (asSlice(u32, interaction.other_end_atoms)) |index| {
        if (!isIndexPresent(index) or index >= atom_count) return .invalid_atom_index;
    }
    return null;
}

/// Input pointers are borrowed only for this call; see coordgen_abi.h. Every
/// validation error api.Input.validate() can produce is reachable here
/// through the C representation. `too_many_items` is not: every span's `len`
/// is already a `uint32_t` in the ABI, so the u32-overflow case
/// api.Input.validate() guards against for the wider Zig `usize` slices
/// cannot occur through this entry point. The discrete/continuous
/// coordinate-generation pipeline (the topology, layout, optimize, and
/// generator layers) is not wired yet, so a structurally valid input
/// deliberately returns `unsupported` rather than a silently-empty result;
/// see the epic's generator beads under cgz-7v2.4.
export fn coordgen_generate(input: *const Input, result: *Result) callconv(.c) core.errors.ErrorCode {
    result.* = .{};

    const atom_count: usize = input.atoms.len;
    if (atom_count == 0) return .empty_graph;

    if (!std.math.isFinite(input.options.precision) or input.options.precision <= 0) {
        return .invalid_option;
    }

    for (asSlice(AtomInput, input.atoms)) |atom| {
        if (validateAtomField(atom, atom_count)) |err| return err;
    }
    for (asSlice(BondInput, input.bonds)) |bond| {
        if (validateBondField(bond, atom_count)) |err| return err;
    }
    for (asSlice(BondInput, input.extra_bonds)) |bond| {
        if (validateBondField(bond, atom_count)) |err| return err;
    }
    for (asSlice(ResidueInput, input.residues)) |residue| {
        if (validateResidueField(residue, atom_count)) |err| return err;
    }
    for (asSlice(ResidueInteractionInput, input.residue_interactions)) |interaction| {
        if (validateResidueInteractionField(interaction, atom_count)) |err| return err;
    }

    return .unsupported;
}

const ResultOwner = struct {
    allocator: std.mem.Allocator,
};

/// A zeroed result (every span null, owner null - the state coordgen_generate
/// leaves on any failure) is always a safe no-op to free, matching
/// coordgen_abi.h's documented "on failure, result is zeroed and requires no
/// cleanup" contract: freeing it anyway must not double-free or crash.
export fn coordgen_result_free(result: *Result) callconv(.c) void {
    if (result.owner) |owner_ptr| {
        const owner: *ResultOwner = @ptrCast(@alignCast(owner_ptr));
        const allocator = owner.allocator;
        if (result.coordinates.ptr) |p| allocator.free(p[0..result.coordinates.len]);
        if (result.input_to_internal.ptr) |p| allocator.free(p[0..result.input_to_internal.len]);
        if (result.internal_to_input.ptr) |p| allocator.free(p[0..result.internal_to_input.len]);
        if (result.effective_bond_orders.ptr) |p| allocator.free(p[0..result.effective_bond_orders.len]);
        if (result.bond_displays.ptr) |p| allocator.free(p[0..result.bond_displays.len]);
        if (result.atom_stereo.ptr) |p| allocator.free(p[0..result.atom_stereo.len]);
        allocator.destroy(owner);
    }
    result.* = .{};
}
test "coordgen_generate rejects an empty graph and leaves result zeroed" {
    const input: Input = .{};
    var result: Result = undefined;
    try std.testing.expectEqual(core.errors.ErrorCode.empty_graph, coordgen_generate(&input, &result));
    try std.testing.expect(result.owner == null);
    try std.testing.expect(result.coordinates.ptr == null);
}

test "coordgen_generate reaches every documented error path that api.Input.validate() can produce" {
    const Case = struct {
        atoms: []const AtomInput,
        bonds: []const BondInput = &.{},
        residues: []const ResidueInput = &.{},
        residue_interactions: []const ResidueInteractionInput = &.{},
        precision: f32 = 1,
        expect: core.errors.ErrorCode,
    };
    const cases = [_]Case{
        .{
            .atoms = &.{.{ .atomic_number = 0 }},
            .expect = .invalid_atomic_number,
        },
        .{
            .atoms = &.{.{ .atomic_number = 119 }},
            .expect = .invalid_atomic_number,
        },
        .{
            .atoms = &.{ .{}, .{} },
            .bonds = &.{.{ .start = 0, .end = 0, .order = 1 }},
            .expect = .invalid_atom_index,
        },
        .{
            .atoms = &.{ .{}, .{} },
            .bonds = &.{.{ .start = 0, .end = 5, .order = 1 }},
            .expect = .invalid_atom_index,
        },
        .{
            .atoms = &.{ .{}, .{} },
            .bonds = &.{.{ .start = 0, .end = 1, .order = 7 }},
            .expect = .invalid_bond_order,
        },
        .{
            .atoms = &.{.{ .stereo = 9 }},
            .expect = .invalid_stereo,
        },
        .{
            // clockwise requires all three references present.
            .atoms = &.{.{ .stereo = 1 }},
            .expect = .invalid_stereo,
        },
        .{
            .atoms = &.{.{ .flags = 1 << 31 }},
            .expect = .invalid_option,
        },
        .{
            .atoms = &.{.{
                .flags = AtomFlags.has_template_coordinates,
                .template_coordinates = .{ .x = std.math.nan(f32), .y = 0 },
            }},
            .expect = .invalid_coordinate,
        },
        .{
            .atoms = &.{.{}},
            .precision = 0,
            .expect = .invalid_option,
        },
        .{
            .atoms = &.{ .{}, .{} },
            .residues = &.{.{ .atom = 9, .closest_ligand_atom = core.ids.invalid_input_index }},
            .expect = .invalid_atom_index,
        },
        .{
            .atoms = &.{ .{}, .{} },
            .residue_interactions = &.{.{ .start = 0, .end = 9 }},
            .expect = .invalid_atom_index,
        },
        .{
            // Structurally sound input: no generator is wired yet.
            .atoms = &.{ .{}, .{} },
            .bonds = &.{.{ .start = 0, .end = 1, .order = 1 }},
            .expect = .unsupported,
        },
    };

    for (cases) |case| {
        const input: Input = .{
            .options = .{ .precision = case.precision },
            .atoms = .{ .ptr = case.atoms.ptr, .len = @intCast(case.atoms.len) },
            .bonds = .{ .ptr = case.bonds.ptr, .len = @intCast(case.bonds.len) },
            .residues = .{ .ptr = case.residues.ptr, .len = @intCast(case.residues.len) },
            .residue_interactions = .{
                .ptr = case.residue_interactions.ptr,
                .len = @intCast(case.residue_interactions.len),
            },
        };
        var result: Result = undefined;
        const code = coordgen_generate(&input, &result);
        try std.testing.expectEqual(case.expect, code);
        if (code != .unsupported) {
            try std.testing.expect(result.owner == null);
        }
    }
}

test "coordgen_result_free is a safe no-op on a zeroed (failure) result" {
    var result: Result = .{};
    coordgen_result_free(&result);
    try std.testing.expect(result.owner == null);
}

test "coordgen_result_free releases every owned span through the stored allocator" {
    const allocator = std.testing.allocator;
    const coordinates = try allocator.alloc(core.math.Vec2, 2);
    const atom_stereo = try allocator.alloc(u32, 2);
    const owner = try allocator.create(ResultOwner);
    owner.* = .{ .allocator = allocator };

    var result: Result = .{
        .coordinates = .{ .ptr = coordinates.ptr, .len = @intCast(coordinates.len) },
        .atom_stereo = .{ .ptr = atom_stereo.ptr, .len = @intCast(atom_stereo.len) },
        .owner = owner,
    };
    coordgen_result_free(&result);
    // std.testing.allocator fails the test binary on any leak at process
    // exit; reaching this line with every span actually freed above (not
    // merely zeroed) is what that leak check is verifying.
    try std.testing.expect(result.owner == null);
    try std.testing.expect(result.coordinates.ptr == null);
}
