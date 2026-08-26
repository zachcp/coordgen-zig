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
const api = @import("api");

/// The installed library is std-only by build policy, so the C entry points
/// cannot reach for libc's allocator. smp_allocator is thread-safe, which the
/// ABI requires: coordgen_abi.h documents concurrent generation on distinct
/// inputs, and the allocator is the only state these entry points share.
const generation_allocator = std.heap.smp_allocator;

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
    // The oracle adapter rejects an out-of-range display with
    // COORDGEN_ERROR_INVALID_STEREO (conformance/oracle_adapter.cpp:206), and
    // a consumer of the shared header must not be able to tell the two
    // implementations apart by which code comes back.
    if (core.chemistry.BondDisplay.fromPublic(bond.display) == null) return .invalid_stereo;
    if (validateBondStereoField(bond, atom_count)) |err| return err;
    return null;
}

fn validFlag(value: u32) bool {
    return value <= 1;
}

/// Mirrors the oracle adapter's option contract exactly
/// (conformance/oracle_adapter.cpp:146-156): a boolean option carries 0 or 1
/// and nothing else, so a caller passing 2 learns it now rather than getting
/// whichever meaning the implementation happened to infer.
fn validateOptions(options: c_abi.Options) ?core.errors.ErrorCode {
    if (!std.math.isFinite(options.precision) or options.precision <= 0) return .invalid_option;
    if (!validFlag(options.score_residue_interactions) or
        !validFlag(options.treat_nonterminal_bonds_to_metal_as_zero_order) or
        !validFlag(options.even_angles) or
        !validFlag(options.skip_minimization) or
        !validFlag(options.force_open_macrocycles) or
        !validFlag(options.constrain_all_atoms) or
        !validFlag(options.reserved) or
        !validFlag(options.debug_coordinates) or
        !validFlag(options.load_templates))
    {
        return .invalid_option;
    }
    // The reserved slot has no safe-API counterpart to carry it: cgz-r25
    // dropped the field from api.Options rather than keep a name for a value
    // that can only be zero. So the C layer is the only place this can be
    // rejected, and it must reject it exactly as the oracle adapter does
    // (conformance/oracle_adapter.cpp:164) - invalid_option for a value that
    // is not a flag at all, unsupported for a well-formed 1.
    if (options.reserved != 0) return .unsupported;
    return null;
}

fn optionalIndex(value: u32) ?core.ids.InputIndex {
    return if (isIndexPresent(value)) value else null;
}

fn stringSlice(view: c_abi.StringView) []const u8 {
    const len: usize = @intCast(view.len);
    return if (view.ptr) |p| p[0..len] else &[_]u8{};
}

/// Borrowed C spans become borrowed-plus-owned safe values: the arrays below
/// are allocated for the duration of the call, while every string and index
/// span inside them still points at caller memory, which the safe API also
/// only borrows.
fn convertAtoms(allocator: std.mem.Allocator, span: []const AtomInput) core.errors.Error![]api.AtomInput {
    const converted = allocator.alloc(api.AtomInput, span.len) catch return error.OutOfMemory;
    errdefer allocator.free(converted);
    for (span, converted) |atom, *destination| {
        destination.* = .{
            .atomic_number = core.chemistry.AtomicNumber.fromPublic(atom.atomic_number) orelse {
                return error.InvalidAtomicNumber;
            },
            .formal_charge = atom.formal_charge,
            .fixed = atom.flags & AtomFlags.fixed != 0,
            .constrained = atom.flags & AtomFlags.constrained != 0,
            .hidden = atom.flags & AtomFlags.hidden != 0,
            .template_coordinates = if (atom.flags & AtomFlags.has_template_coordinates != 0)
                atom.template_coordinates
            else
                null,
            .coordinates_3d = if (atom.flags & AtomFlags.has_coordinates_3d != 0)
                atom.coordinates_3d
            else
                null,
            .stereo = .{
                .value = core.chemistry.AtomStereo.fromPublic(atom.stereo) orelse return error.InvalidStereo,
                .looking_from = optionalIndex(atom.stereo_looking_from),
                .atom_a = optionalIndex(atom.stereo_atom_a),
                .atom_b = optionalIndex(atom.stereo_atom_b),
            },
        };
    }
    return converted;
}

fn convertBonds(allocator: std.mem.Allocator, span: []const BondInput) core.errors.Error![]api.BondInput {
    const converted = allocator.alloc(api.BondInput, span.len) catch return error.OutOfMemory;
    errdefer allocator.free(converted);
    for (span, converted) |bond, *destination| {
        destination.* = .{
            .start = bond.start,
            .end = bond.end,
            .order = core.chemistry.BondOrder.fromInt(bond.order) orelse return error.InvalidBondOrder,
            .skip = bond.flags & BondFlags.skip != 0,
            .stereo = .{
                .value = core.chemistry.BondStereo.fromPublic(bond.stereo) orelse return error.InvalidStereo,
                .atom_a = optionalIndex(bond.stereo_atom_a),
                .atom_b = optionalIndex(bond.stereo_atom_b),
            },
            .display = core.chemistry.BondDisplay.fromPublic(bond.display) orelse {
                return error.InvalidStereo;
            },
            .crossing_penalty_multiplier = bond.crossing_penalty_multiplier,
        };
    }
    return converted;
}

fn convertResidues(allocator: std.mem.Allocator, span: []const ResidueInput) core.errors.Error![]api.ResidueInput {
    const converted = allocator.alloc(api.ResidueInput, span.len) catch return error.OutOfMemory;
    errdefer allocator.free(converted);
    for (span, converted) |residue, *destination| {
        destination.* = .{
            .atom = residue.atom,
            .chain = stringSlice(residue.chain),
            .residue_number = residue.residue_number,
            .closest_ligand_atom = optionalIndex(residue.closest_ligand_atom),
        };
    }
    return converted;
}

fn convertResidueInteractions(
    allocator: std.mem.Allocator,
    span: []const ResidueInteractionInput,
) core.errors.Error![]api.ResidueInteractionInput {
    const converted = allocator.alloc(api.ResidueInteractionInput, span.len) catch return error.OutOfMemory;
    errdefer allocator.free(converted);
    for (span, converted) |interaction, *destination| {
        destination.* = .{
            .start = interaction.start,
            .end = interaction.end,
            .other_start_atoms = asSlice(u32, interaction.other_start_atoms),
            .other_end_atoms = asSlice(u32, interaction.other_end_atoms),
            .crossing_penalty_multiplier = interaction.crossing_penalty_multiplier,
        };
    }
    return converted;
}

fn encodeEnumSpan(allocator: std.mem.Allocator, values: anytype) core.errors.Error![]u32 {
    const encoded = allocator.alloc(u32, values.len) catch return error.OutOfMemory;
    for (values, encoded) |value, *destination| destination.* = @backingInt(value);
    return encoded;
}

fn vec2Span(slice: []core.math.Vec2) c_abi.Vec2Span {
    return .{ .ptr = slice.ptr, .len = @intCast(slice.len) };
}

fn u32Span(slice: []u32) c_abi.U32Span {
    return .{ .ptr = slice.ptr, .len = @intCast(slice.len) };
}

/// Moves `generated` into the C result. The three enum-typed observables are
/// re-encoded as u32 spans because the ABI carries them that way; the other
/// three slices are handed over unchanged. On success `generated` no longer
/// owns those three, which is why the caller releases exactly the re-encoded
/// originals and never calls `api.Result.deinit`.
fn publishResult(
    allocator: std.mem.Allocator,
    generated: api.Result,
    result: *Result,
) core.errors.Error!void {
    const effective_bond_orders = try encodeEnumSpan(allocator, generated.effective_bond_orders);
    errdefer allocator.free(effective_bond_orders);
    const bond_displays = try encodeEnumSpan(allocator, generated.bond_displays);
    errdefer allocator.free(bond_displays);
    const atom_stereo = try encodeEnumSpan(allocator, generated.atom_stereo);
    errdefer allocator.free(atom_stereo);
    const owner = allocator.create(ResultOwner) catch return error.OutOfMemory;
    owner.* = .{ .allocator = allocator };

    // Nothing below can fail, so the ownership move is atomic from the
    // caller's side: either every span is published, or none is.
    result.* = .{
        .coordinates = vec2Span(generated.coordinates),
        .input_to_internal = u32Span(generated.input_to_internal),
        .internal_to_input = u32Span(generated.internal_to_input),
        .effective_bond_orders = u32Span(effective_bond_orders),
        .bond_displays = u32Span(bond_displays),
        .atom_stereo = u32Span(atom_stereo),
        .clean_pose = @intFromBool(generated.clean_pose),
        .owner = owner,
    };
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
/// cannot occur through this entry point.
///
/// Coverage is partial by domain, never by size, and a domain no native phase
/// owns yet returns `unsupported` rather than a silently-empty success. The
/// rejected set is enumerated on `api.generate`.
export fn coordgen_generate(input: *const Input, result: *Result) callconv(.c) core.errors.ErrorCode {
    result.* = .{};

    const atom_count: usize = input.atoms.len;
    if (atom_count == 0) return .empty_graph;

    if (validateOptions(input.options)) |err| return err;

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

    return generateThroughSafeApi(generation_allocator, input, result);
}

fn generateThroughSafeApi(
    allocator: std.mem.Allocator,
    input: *const Input,
    result: *Result,
) core.errors.ErrorCode {
    const atoms = convertAtoms(allocator, asSlice(AtomInput, input.atoms)) catch |err| {
        return core.errors.code(err);
    };
    defer allocator.free(atoms);
    const bonds = convertBonds(allocator, asSlice(BondInput, input.bonds)) catch |err| {
        return core.errors.code(err);
    };
    defer allocator.free(bonds);
    const extra_bonds = convertBonds(allocator, asSlice(BondInput, input.extra_bonds)) catch |err| {
        return core.errors.code(err);
    };
    defer allocator.free(extra_bonds);
    const residues = convertResidues(allocator, asSlice(ResidueInput, input.residues)) catch |err| {
        return core.errors.code(err);
    };
    defer allocator.free(residues);
    const residue_interactions = convertResidueInteractions(
        allocator,
        asSlice(ResidueInteractionInput, input.residue_interactions),
    ) catch |err| return core.errors.code(err);
    defer allocator.free(residue_interactions);

    const template_directory = input.options.template_directory;
    const safe_input = api.Input{
        .atoms = atoms,
        .bonds = bonds,
        .residues = residues,
        .residue_interactions = residue_interactions,
        .extra_bonds = extra_bonds,
        .options = .{
            .precision = input.options.precision,
            .score_residue_interactions = input.options.score_residue_interactions != 0,
            .treat_nonterminal_bonds_to_metal_as_zero_order = input.options.treat_nonterminal_bonds_to_metal_as_zero_order != 0,
            .even_angles = input.options.even_angles != 0,
            .skip_minimization = input.options.skip_minimization != 0,
            .force_open_macrocycles = input.options.force_open_macrocycles != 0,
            .constrain_all_atoms = input.options.constrain_all_atoms != 0,
            .debug_coordinates = input.options.debug_coordinates != 0,
            .load_templates = input.options.load_templates != 0,
            // A present-but-empty view is still a request for process-global
            // template state, which no phase owns; only an absent view selects
            // the built-in data. The oracle adapter refuses the same shape.
            .template_directory = if (template_directory.ptr != null or template_directory.len != 0)
                stringSlice(template_directory)
            else
                null,
        },
    };

    var generated = api.generate(allocator, safe_input) catch |err| return core.errors.code(err);
    publishResult(allocator, generated, result) catch |err| {
        generated.deinit();
        result.* = .{};
        return core.errors.code(err);
    };
    // publishResult moved the other three slices into `result`; these are the
    // originals of the three it re-encoded.
    allocator.free(generated.atom_stereo);
    allocator.free(generated.bond_displays);
    allocator.free(generated.effective_bond_orders);
    return .ok;
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
        even_angles: u32 = 0,
        reserved: u32 = 0,
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
            // A boolean option carrying a value other than 0 or 1.
            .atoms = &.{.{}},
            .expect = .invalid_option,
            .even_angles = 2,
        },
        .{
            .atoms = &.{ .{}, .{} },
            .bonds = &.{.{ .start = 0, .end = 1, .order = 1, .display = 9 }},
            .expect = .invalid_stereo,
        },
        .{
            // A domain no native phase owns: rejected by name, never as an
            // empty success. Hidden atoms are excluded from canonical order
            // upstream, and no native phase reproduces that yet.
            .atoms = &.{ .{}, .{ .flags = AtomFlags.hidden } },
            .bonds = &.{.{ .start = 0, .end = 1, .order = 1 }},
            .expect = .unsupported,
        },
        .{
            // cgz-r25's reserved slot. A non-flag value is invalid_option;
            // a well-formed 1 is unsupported, matching the oracle adapter.
            // api.Options cannot carry this, so this boundary is the only
            // place it can be caught.
            .atoms = &.{.{}},
            .reserved = 2,
            .expect = .invalid_option,
        },
        .{
            .atoms = &.{ .{}, .{} },
            .bonds = &.{.{ .start = 0, .end = 1, .order = 1 }},
            .reserved = 1,
            .expect = .unsupported,
        },
        .{
            // Structurally sound and in-domain.
            .atoms = &.{ .{}, .{} },
            .bonds = &.{.{ .start = 0, .end = 1, .order = 1 }},
            .expect = .ok,
        },
    };

    for (cases) |case| {
        const input: Input = .{
            .options = .{
                .precision = case.precision,
                .even_angles = case.even_angles,
                .reserved = case.reserved,
            },
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
        if (code == .ok) {
            try std.testing.expect(result.owner != null);
            coordgen_result_free(&result);
        } else {
            // Every failure path leaves the documented zero value, so a caller
            // that frees unconditionally is still correct.
            try std.testing.expect(result.owner == null);
            try std.testing.expect(result.coordinates.ptr == null);
        }
    }
}

test "coordgen_generate publishes every observable span and frees them as one unit" {
    const atoms = [_]AtomInput{ .{}, .{}, .{} };
    const bonds = [_]BondInput{
        .{ .start = 0, .end = 1, .order = 1 },
        .{ .start = 1, .end = 2, .order = 2 },
    };
    const input: Input = .{
        .atoms = .{ .ptr = &atoms, .len = atoms.len },
        .bonds = .{ .ptr = &bonds, .len = bonds.len },
    };
    var result: Result = undefined;
    try std.testing.expectEqual(core.errors.ErrorCode.ok, coordgen_generate(&input, &result));

    try std.testing.expectEqual(@as(u32, atoms.len), result.coordinates.len);
    try std.testing.expectEqual(@as(u32, atoms.len), result.input_to_internal.len);
    try std.testing.expectEqual(@as(u32, atoms.len), result.internal_to_input.len);
    try std.testing.expectEqual(@as(u32, atoms.len), result.atom_stereo.len);
    try std.testing.expectEqual(@as(u32, bonds.len), result.effective_bond_orders.len);
    try std.testing.expectEqual(@as(u32, bonds.len), result.bond_displays.len);

    const coordinates = result.coordinates.ptr.?[0..result.coordinates.len];
    for (coordinates) |position| try std.testing.expect(position.isFinite());
    for (bonds) |bond| {
        const delta_x = coordinates[bond.start].x - coordinates[bond.end].x;
        const delta_y = coordinates[bond.start].y - coordinates[bond.end].y;
        try std.testing.expectApproxEqAbs(
            c_abi.bond_length,
            @sqrt(delta_x * delta_x + delta_y * delta_y),
            0.001,
        );
    }
    // The enum observables cross the ABI as their pinned numeric values.
    const effective = result.effective_bond_orders.ptr.?[0..result.effective_bond_orders.len];
    try std.testing.expectEqualSlices(u32, &.{ 1, 2 }, effective);
    const displays = result.bond_displays.ptr.?[0..result.bond_displays.len];
    try std.testing.expectEqualSlices(u32, &.{ 0, 0 }, displays);

    coordgen_result_free(&result);
    try std.testing.expect(result.owner == null);
    try std.testing.expect(result.coordinates.ptr == null);
    try std.testing.expect(result.atom_stereo.ptr == null);
}

test "coordgen_generate is repeatable and independent across concurrent contexts" {
    const atoms = [_]AtomInput{ .{}, .{}, .{}, .{} };
    const bonds = [_]BondInput{
        .{ .start = 0, .end = 1, .order = 1 },
        .{ .start = 1, .end = 2, .order = 1 },
        .{ .start = 2, .end = 3, .order = 1 },
    };
    const input: Input = .{
        .atoms = .{ .ptr = &atoms, .len = atoms.len },
        .bonds = .{ .ptr = &bonds, .len = bonds.len },
    };

    const Worker = struct {
        fn run(shared: *const Input, out: *[4]core.math.Vec2, ok: *bool) void {
            var local: Result = undefined;
            ok.* = coordgen_generate(shared, &local) == .ok;
            if (!ok.*) return;
            @memcpy(out, local.coordinates.ptr.?[0..4]);
            coordgen_result_free(&local);
        }
    };

    var first: [4]core.math.Vec2 = undefined;
    var second: [4]core.math.Vec2 = undefined;
    var first_ok = false;
    var second_ok = false;
    const left = try std.Thread.spawn(.{}, Worker.run, .{ &input, &first, &first_ok });
    const right = try std.Thread.spawn(.{}, Worker.run, .{ &input, &second, &second_ok });
    left.join();
    right.join();

    try std.testing.expect(first_ok);
    try std.testing.expect(second_ok);
    try std.testing.expectEqualSlices(core.math.Vec2, &first, &second);
}

/// The seam `coordgen_generate` delegates to, driven with an injectable
/// allocator. The export itself pins `generation_allocator` because the
/// installed library is std-only and the ABI promises thread safety, so it
/// cannot take one from the caller; everything it does before delegating is
/// field validation that allocates nothing. Sweeping the seam is therefore
/// sweeping the entry point.
fn generateThroughSafeApiAndDiscard(allocator: std.mem.Allocator, input: *const Input) !void {
    var result: Result = .{};
    switch (generateThroughSafeApi(allocator, input, &result)) {
        .ok => {},
        .out_of_memory => {
            // The contract in coordgen_abi.h is that a failed generation
            // leaves the result zeroed and needing no cleanup. Assert that
            // rather than assuming it, since a partially published result
            // escaping here would leak through the C boundary where no
            // errdefer can reach it.
            if (result.owner != null) return error.PartialResultEscaped;
            return error.OutOfMemory;
        },
        else => |code| {
            std.debug.print("unexpected error code: {t}\n", .{code});
            return error.UnexpectedErrorCode;
        },
    }
    coordgen_result_free(&result);
}

test "coordgen_generate reports and cleans up failure at every allocation index" {
    const atoms = [_]AtomInput{ .{}, .{}, .{}, .{} };
    // Atom 3 carries the residue and must therefore have no structural bond:
    // a residue representative is placed by the residue phase, not laid out.
    const bonds = [_]BondInput{
        .{ .start = 0, .end = 1, .order = 1 },
        .{ .start = 1, .end = 2, .order = 2 },
    };
    const residues = [_]ResidueInput{
        .{ .atom = 3, .chain = .{ .ptr = "A", .len = 1 }, .residue_number = 9, .closest_ligand_atom = 1 },
    };
    const interactions = [_]ResidueInteractionInput{.{ .start = 3, .end = 1 }};
    // Residues and interactions are included so the sweep reaches the
    // conversion allocations that a bare atoms-and-bonds input never makes.
    const input: Input = .{
        .atoms = .{ .ptr = &atoms, .len = atoms.len },
        .bonds = .{ .ptr = &bonds, .len = bonds.len },
        .residues = .{ .ptr = &residues, .len = residues.len },
        .residue_interactions = .{ .ptr = &interactions, .len = interactions.len },
    };

    // Discharge the first-call warm-up, so the counting run inside
    // checkAllocationFailures does not overcount and leave the last index
    // unexercised.
    try generateThroughSafeApiAndDiscard(std.testing.allocator, &input);
    var leak_check = std.testing.FailingAllocator.init(std.testing.allocator, .{});
    try generateThroughSafeApiAndDiscard(leak_check.allocator(), &input);
    try std.testing.expectEqual(leak_check.allocated_bytes, leak_check.freed_bytes);

    try core.oom.checkAllocationFailures(
        std.testing.allocator,
        generateThroughSafeApiAndDiscard,
        .{&input},
    );
}

/// The C ABI's own allocation-site floor. coordgen_generate pins
/// generation_allocator, so the seam it delegates to is the injectable point;
/// everything before that delegation is non-allocating field validation.
const allocation_site_table = @embedFile("allocation_site_floors");

test "the C entry point does not silently stop reaching allocation sites" {
    const atoms = [_]AtomInput{ .{}, .{}, .{}, .{} };
    const bonds = [_]BondInput{
        .{ .start = 0, .end = 1, .order = 1 },
        .{ .start = 1, .end = 2, .order = 2 },
    };
    const residues = [_]ResidueInput{
        .{ .atom = 3, .chain = .{ .ptr = "A", .len = 1 }, .residue_number = 9, .closest_ligand_atom = 1 },
    };
    const interactions = [_]ResidueInteractionInput{.{ .start = 3, .end = 1 }};
    const input: Input = .{
        .atoms = .{ .ptr = &atoms, .len = atoms.len },
        .bonds = .{ .ptr = &bonds, .len = bonds.len },
        .residues = .{ .ptr = &residues, .len = residues.len },
        .residue_interactions = .{ .ptr = &interactions, .len = interactions.len },
    };
    const measured = try core.oom.countAllocationSites(
        std.testing.allocator,
        generateThroughSafeApiAndDiscard,
        .{&input},
    );
    try core.oom.expectSiteFloor(allocation_site_table, "coordgen_generate", "residue ligand", measured);
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

// ---------------------------------------------------------------------------
// Fuzz target (cgz-7v2.4.4)
//
// This lives here rather than in tests/fuzz_targets.zig because the exports
// are `export fn`, not `pub`, so this file is the only place that can call
// them without linking the library and re-declaring the symbols.
//
// It is a SEPARATE surface from the native targets, not a thin wrapper over
// them, and the reason is concrete: every enum on this boundary crosses as a
// bare `u32`. An out-of-range `atomic_number` or `order` is unrepresentable in
// the Zig API - the type system refuses it - and is trivially representable
// here. Those values reach the validators only through this door.
// ---------------------------------------------------------------------------

const fuzz_max_atoms = 10;
const fuzz_max_bonds = 14;

const FuzzInput = struct {
    atoms: [fuzz_max_atoms]AtomInput = undefined,
    bonds: [fuzz_max_bonds]BondInput = undefined,
    atom_count: u32 = 0,
    bond_count: u32 = 0,
    /// The structurally valid lane must reach generation rather than stopping
    /// at a C DTO validator. Generation may still return `unsupported` for a
    /// valid stereo mode the native port does not own yet.
    expect_generation: bool = false,

    fn view(self: *const FuzzInput) Input {
        return .{
            .atoms = .{ .ptr = &self.atoms, .len = self.atom_count },
            .bonds = .{ .ptr = &self.bonds, .len = self.bond_count },
        };
    }
};

/// Each field is drawn mostly in range, occasionally not.
///
/// The first version of this drew `atomic_number` from the whole `u16` and
/// endpoints from the whole `u8`, so nearly every input died at the first
/// validator: 619 runs, 19 unique, 0.29% coverage - a target that tests the
/// validator and nothing behind it. With the weights below, a comparable
/// session reached 650 runs, 43 unique, 0.51%. Out-of-range values have to be
/// the MINORITY for the generator itself to be reached, while still appearing
/// often enough that every rejection path is exercised; these weights are that
/// trade, not a preference for valid input.
///
/// Neither number is evidence of good coverage - both are short sessions, and
/// coverage accumulates across runs. What they show is the direction the
/// rebalance moved, which is the only claim being made.
fn buildFuzzInput(smith: *std.testing.Smith) FuzzInput {
    var built: FuzzInput = .{};
    // Per-field validity compounds into almost no valid whole inputs. Choose
    // the lane once instead: three quarters are connected, validator-clean
    // molecules that must reach generation; one quarter explores raw C values.
    built.expect_generation = smith.boolWeighted(3, 1);
    built.atom_count = if (built.expect_generation)
        smith.valueRangeAtMost(u8, 1, fuzz_max_atoms)
    else
        smith.valueRangeAtMost(u8, 0, fuzz_max_atoms);
    for (built.atoms[0..built.atom_count]) |*atom| {
        atom.* = .{
            // Mostly ordinary elements; sometimes a value above 118, which is
            // representable only here and whose documented answer is
            // `invalid_atomic_number`.
            .atomic_number = if (built.expect_generation or smith.boolWeighted(9, 1))
                smith.valueRangeAtMost(u8, 1, 30)
            else
                smith.value(u16),
            .formal_charge = smith.valueRangeAtMost(i8, -3, 3),
            .flags = if (built.expect_generation) 0 else smith.valueRangeAtMost(u8, 0, 63),
        };
        if (built.expect_generation) {
            atom.stereo = smith.valueRangeAtMost(u8, 0, 4);
            if (atom.stereo == 1 or atom.stereo == 2) {
                atom.stereo_looking_from = @intCast(smith.index(built.atom_count));
                atom.stereo_atom_a = @intCast(smith.index(built.atom_count));
                atom.stereo_atom_b = @intCast(smith.index(built.atom_count));
            }
        } else {
            atom.stereo = smith.valueRangeAtMost(u8, 0, 8);
            atom.stereo_looking_from = smith.value(u8);
            atom.stereo_atom_a = smith.value(u8);
            atom.stereo_atom_b = smith.value(u8);
        }
    }
    if (built.expect_generation) {
        // A connected chain keeps this lane inside the generator's owned
        // domain. Atom fields and stereo still vary on every iteration.
        for (0..built.atom_count - 1) |index| {
            built.bonds[index] = .{
                .start = @intCast(index),
                .end = @intCast(index + 1),
                .order = smith.valueRangeAtMost(u8, 1, 3),
            };
        }
        built.bond_count = built.atom_count - 1;
        return built;
    }
    var count: u32 = 0;
    while (count < fuzz_max_bonds and !smith.eosWeightedSimple(6, 1)) {
        built.bonds[count] = .{
            .start = if (built.atom_count != 0 and smith.boolWeighted(6, 1))
                @intCast(smith.index(built.atom_count))
            else
                smith.value(u8),
            .end = if (built.atom_count != 0 and smith.boolWeighted(6, 1))
                @intCast(smith.index(built.atom_count))
            else
                smith.value(u8),
            // Above 3 is not a bond order at all.
            .order = if (smith.boolWeighted(6, 1))
                smith.valueRangeAtMost(u8, 0, 3)
            else
                smith.valueRangeAtMost(u8, 0, 6),
            .flags = smith.valueRangeAtMost(u8, 0, 7),
            .stereo = smith.valueRangeAtMost(u8, 0, 6),
            .display = smith.valueRangeAtMost(u8, 0, 6),
            .crossing_penalty_multiplier = 1,
        };
        count += 1;
    }
    built.bond_count = count;
    return built;
}

fn cAbiContract(_: void, smith: *std.testing.Smith) anyerror!void {
    var built = buildFuzzInput(smith);
    const input = built.view();
    var result: Result = undefined;
    const code = coordgen_generate(&input, &result);
    if (built.expect_generation) {
        try std.testing.expect(code == .ok or code == .unsupported);
    }

    if (code != .ok) {
        // coordgen_abi.h promises a zeroed result on failure, needing no
        // cleanup. A partially published result escaping here would leak
        // across the C boundary where no errdefer can reach it, and the
        // caller - following the header - would never free it.
        try std.testing.expect(result.owner == null);
        try std.testing.expect(result.coordinates.ptr == null);
        try std.testing.expect(result.input_to_internal.ptr == null);
        try std.testing.expect(result.internal_to_input.ptr == null);
        try std.testing.expect(result.effective_bond_orders.ptr == null);
        try std.testing.expect(result.bond_displays.ptr == null);
        try std.testing.expect(result.atom_stereo.ptr == null);
        // Freeing a zeroed result is documented as a safe no-op, so a caller
        // that frees unconditionally must not be punished for it.
        coordgen_result_free(&result);
        return;
    }

    try std.testing.expectEqual(built.atom_count, result.coordinates.len);
    try std.testing.expectEqual(built.atom_count, result.input_to_internal.len);
    try std.testing.expectEqual(built.atom_count, result.internal_to_input.len);
    try std.testing.expectEqual(built.atom_count, result.atom_stereo.len);
    try std.testing.expectEqual(built.bond_count, result.effective_bond_orders.len);
    try std.testing.expectEqual(built.bond_count, result.bond_displays.len);
    if (built.atom_count != 0) try std.testing.expect(result.owner != null);

    const coordinates = if (result.coordinates.ptr) |ptr|
        ptr[0..result.coordinates.len]
    else
        &[_]core.math.Vec2{};
    for (coordinates) |coordinate| try std.testing.expect(coordinate.isFinite());

    // clean_pose crosses as u32 and is a boolean; anything else means the
    // encoding drifted.
    try std.testing.expect(result.clean_pose == 0 or result.clean_pose == 1);

    coordgen_result_free(&result);
    try std.testing.expect(result.owner == null);
    // Idempotent by contract: the second free sees the zeroed value the first
    // left, which is the same state a failed generation leaves.
    coordgen_result_free(&result);
}

const c_abi_fuzz_seeds = @import("fuzz_seeds").c_abi_contract_seeds;

test "fuzz: c abi contract for atom and bond inputs" {
    try std.testing.fuzz({}, cAbiContract, .{ .corpus = c_abi_fuzz_seeds });
}
