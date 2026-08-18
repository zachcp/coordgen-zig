const std = @import("std");
const core = @import("core");
const model = @import("model");
const canonical = @import("canonical.zig");

pub const Prepared = struct {
    working: model.WorkingGraph,

    pub fn deinit(self: *Prepared) void {
        self.working.deinit();
        self.* = undefined;
    }
};

/// Copies a validated safe-API-shaped input into canonical allocator-owned
/// storage. `anytype` keeps the topology layer below the public API layer while
/// preserving one DTO source of truth; generator passes `api.Input` directly.
pub fn init(allocator: std.mem.Allocator, input: anytype) core.errors.Error!Prepared {
    if (input.atoms.len == 0) return error.EmptyGraph;
    if (input.atoms.len > std.math.maxInt(u32) or
        input.bonds.len > std.math.maxInt(u32) or
        input.extra_bonds.len > std.math.maxInt(u32) or
        input.residues.len > std.math.maxInt(u32) or
        input.residue_interactions.len > std.math.maxInt(u32)) return error.TooManyItems;

    const atom_descriptors = allocator.alloc(canonical.Atom, input.atoms.len) catch return error.OutOfMemory;
    defer allocator.free(atom_descriptors);
    for (input.atoms, atom_descriptors) |atom, *descriptor| {
        descriptor.* = .{ .atomic_number = atom.atomic_number, .hidden = atom.hidden };
    }

    const degrees = allocator.alloc(u32, input.atoms.len) catch return error.OutOfMemory;
    defer allocator.free(degrees);
    @memset(degrees, 0);
    for (input.bonds) |bond| {
        try validateEndpoints(bond.start, bond.end, input.atoms.len);
        degrees[bond.start] += 1;
        degrees[bond.end] += 1;
    }
    for (input.extra_bonds) |bond| try validateEndpoints(bond.start, bond.end, input.atoms.len);

    const bond_descriptors = allocator.alloc(canonical.Bond, input.bonds.len) catch return error.OutOfMemory;
    defer allocator.free(bond_descriptors);
    for (input.bonds, bond_descriptors) |bond, *descriptor| {
        descriptor.* = .{
            .start = bond.start,
            .end = bond.end,
            .effective_order = effectiveOrder(
                input.options.treat_nonterminal_bonds_to_metal_as_zero_order,
                atom_descriptors,
                degrees,
                bond,
            ),
            .skip = bond.skip,
        };
    }

    var ordering = try canonical.order(allocator, atom_descriptors, bond_descriptors);
    defer ordering.deinit();
    var order_map = try model.InputOrderMap.init(allocator, ordering.internal_to_input);
    errdefer order_map.deinit();

    const atoms = allocator.alloc(model.Atom, input.atoms.len) catch return error.OutOfMemory;
    errdefer allocator.free(atoms);
    for (ordering.internal_to_input, 0..) |input_index, internal_index| {
        const source = input.atoms[input_index];
        atoms[internal_index] = .{
            .id = core.ids.AtomId.fromIndex(@intCast(internal_index)),
            .input_index = input_index,
            .atomic_number = source.atomic_number,
            .formal_charge = source.formal_charge,
            .fixed = source.fixed,
            .constrained = source.constrained,
            .hidden = source.hidden,
            .template_coordinates = source.template_coordinates,
            .coordinates_3d = source.coordinates_3d,
            .stereo = source.stereo.value,
            .stereo_looking_from = remapOptional(order_map, source.stereo.looking_from),
            .stereo_atom_a = remapOptional(order_map, source.stereo.atom_a),
            .stereo_atom_b = remapOptional(order_map, source.stereo.atom_b),
        };
    }

    const bonds = allocator.alloc(model.Bond, input.bonds.len) catch return error.OutOfMemory;
    errdefer allocator.free(bonds);
    const emitted = allocator.alloc(bool, input.bonds.len) catch return error.OutOfMemory;
    defer allocator.free(emitted);
    @memset(emitted, false);
    var bond_internal_index: usize = 0;
    for (ordering.structural_bonds) |input_index| {
        bonds[bond_internal_index] = makeBond(
            @intCast(bond_internal_index),
            input_index,
            input.bonds[input_index],
            bond_descriptors[input_index].effective_order,
            order_map,
        );
        emitted[input_index] = true;
        bond_internal_index += 1;
    }
    for (input.bonds, 0..) |bond, input_index| {
        if (emitted[input_index]) continue;
        bonds[bond_internal_index] = makeBond(
            @intCast(bond_internal_index),
            @intCast(input_index),
            bond,
            bond_descriptors[input_index].effective_order,
            order_map,
        );
        bond_internal_index += 1;
    }

    const extra_bonds = allocator.alloc(model.ExtraBond, input.extra_bonds.len) catch return error.OutOfMemory;
    errdefer allocator.free(extra_bonds);
    for (input.extra_bonds, extra_bonds, 0..) |source, *destination, input_index| {
        destination.* = .{
            .input_index = @intCast(input_index),
            .start = remap(order_map, source.start),
            .end = remap(order_map, source.end),
            .input_order = source.order,
            .effective_order = source.order,
            .skip = source.skip,
            .stereo = source.stereo.value,
            .stereo_atom_a = remapOptional(order_map, source.stereo.atom_a),
            .stereo_atom_b = remapOptional(order_map, source.stereo.atom_b),
            .display = source.display,
            .crossing_penalty_multiplier = source.crossing_penalty_multiplier,
        };
    }

    var proximity_count = countProximity(bonds) + input.residue_interactions.len;
    for (extra_bonds) |bond| proximity_count += @intFromBool(!bond.skip and bond.effective_order == .zero);
    const proximity_relations = allocator.alloc(model.ProximityRelation, proximity_count) catch return error.OutOfMemory;
    errdefer allocator.free(proximity_relations);
    var proximity_index: usize = 0;
    for (bonds) |bond| {
        if (bond.skip or bond.effective_order != .zero) continue;
        proximity_relations[proximity_index] = .{ .bond = bond.id };
        proximity_index += 1;
    }
    for (extra_bonds) |bond| {
        if (bond.skip or bond.effective_order != .zero) continue;
        proximity_relations[proximity_index] = .{ .extra_bond = bond.input_index };
        proximity_index += 1;
    }
    for (input.residue_interactions, 0..) |_, interaction_index| {
        proximity_relations[proximity_index] = .{ .residue_interaction = @intCast(interaction_index) };
        proximity_index += 1;
    }

    const residues = allocator.alloc(model.Residue, input.residues.len) catch return error.OutOfMemory;
    errdefer allocator.free(residues);
    var string_len: usize = 0;
    for (input.residues) |residue| {
        if (residue.atom >= input.atoms.len or
            (residue.closest_ligand_atom != null and residue.closest_ligand_atom.? >= input.atoms.len))
            return error.InvalidAtomIndex;
        string_len += residue.chain.len;
        if (string_len > std.math.maxInt(u32)) return error.TooManyItems;
    }
    const string_bytes = allocator.alloc(u8, string_len) catch return error.OutOfMemory;
    errdefer allocator.free(string_bytes);
    var string_offset: usize = 0;
    for (input.residues, residues, 0..) |source, *destination, index| {
        @memcpy(string_bytes[string_offset..][0..source.chain.len], source.chain);
        destination.* = .{
            .id = core.ids.ResidueId.fromIndex(@intCast(index)),
            .atom = remap(order_map, source.atom),
            .chain_start = @intCast(string_offset),
            .chain_len = @intCast(source.chain.len),
            .residue_number = source.residue_number,
            .closest_ligand_atom = remapOptional(order_map, source.closest_ligand_atom),
        };
        string_offset += source.chain.len;
    }

    var interaction_atom_count: usize = 0;
    for (input.residue_interactions) |interaction| {
        try validateIndex(interaction.start, input.atoms.len);
        try validateIndex(interaction.end, input.atoms.len);
        interaction_atom_count += interaction.other_start_atoms.len + interaction.other_end_atoms.len;
        if (interaction_atom_count > std.math.maxInt(u32)) return error.TooManyItems;
    }
    const interaction_atoms = allocator.alloc(core.ids.AtomId, interaction_atom_count) catch return error.OutOfMemory;
    errdefer allocator.free(interaction_atoms);
    const interactions = allocator.alloc(model.ResidueInteraction, input.residue_interactions.len) catch return error.OutOfMemory;
    errdefer allocator.free(interactions);
    var interaction_offset: usize = 0;
    for (input.residue_interactions, interactions) |source, *destination| {
        const start_offset = interaction_offset;
        for (source.other_start_atoms) |atom| {
            if (atom >= input.atoms.len) return error.InvalidAtomIndex;
            interaction_atoms[interaction_offset] = remap(order_map, atom);
            interaction_offset += 1;
        }
        const end_offset = interaction_offset;
        for (source.other_end_atoms) |atom| {
            if (atom >= input.atoms.len) return error.InvalidAtomIndex;
            interaction_atoms[interaction_offset] = remap(order_map, atom);
            interaction_offset += 1;
        }
        destination.* = .{
            .start = remap(order_map, source.start),
            .end = remap(order_map, source.end),
            .other_start_atoms = .{ .start = @intCast(start_offset), .len = @intCast(end_offset - start_offset) },
            .other_end_atoms = .{ .start = @intCast(end_offset), .len = @intCast(interaction_offset - end_offset) },
            .crossing_penalty_multiplier = source.crossing_penalty_multiplier,
        };
    }

    try flagCrossAtoms(allocator, atoms, bonds);

    return .{ .working = .{
        .allocator = allocator,
        .atoms = atoms,
        .bonds = bonds,
        .extra_bonds = extra_bonds,
        .residues = residues,
        .residue_interactions = interactions,
        .residue_interaction_atoms = interaction_atoms,
        .string_bytes = string_bytes,
        .order = order_map,
        .active_atom_count = ordering.visible_count,
        .proximity_relations = proximity_relations,
    } };
}

fn effectiveOrder(enabled: bool, atoms: []const canonical.Atom, degrees: []const u32, bond: anytype) core.chemistry.BondOrder {
    if (!enabled or bond.skip or (bond.order != .single and bond.order != .double)) return bond.order;
    if (degrees[bond.start] == 1 or degrees[bond.end] == 1) return bond.order;
    if (isMetal(atoms[bond.start].atomic_number) or isMetal(atoms[bond.end].atomic_number)) return .zero;
    return bond.order;
}

pub fn isMetal(number: core.chemistry.AtomicNumber) bool {
    const value = @backingInt(number);
    return (value >= 3 and value <= 4) or (value >= 11 and value <= 13) or
        (value >= 19 and value <= 32) or (value >= 37 and value <= 51) or
        (value >= 55 and value <= 84) or
        (value >= 87 and value <= 112);
}

fn makeBond(internal_index: u32, input_index: u32, source: anytype, effective_order: core.chemistry.BondOrder, order_map: model.InputOrderMap) model.Bond {
    return .{
        .id = core.ids.BondId.fromIndex(internal_index),
        .input_index = input_index,
        .start = remap(order_map, source.start),
        .end = remap(order_map, source.end),
        .input_order = source.order,
        .effective_order = effective_order,
        .skip = source.skip,
        .stereo = source.stereo.value,
        .stereo_atom_a = remapOptional(order_map, source.stereo.atom_a),
        .stereo_atom_b = remapOptional(order_map, source.stereo.atom_b),
        .display = source.display,
        .crossing_penalty_multiplier = source.crossing_penalty_multiplier,
    };
}

fn remap(order_map: model.InputOrderMap, input_index: core.ids.InputIndex) core.ids.AtomId {
    return core.ids.AtomId.fromIndex(order_map.input_to_internal[input_index]);
}

fn remapOptional(order_map: model.InputOrderMap, input_index: ?core.ids.InputIndex) core.ids.AtomId {
    return if (input_index) |index| remap(order_map, index) else .invalid;
}

fn validateEndpoints(start: u32, end: u32, count: usize) core.errors.Error!void {
    if (start >= count or end >= count or start == end) return error.InvalidAtomIndex;
}

fn validateIndex(index: u32, count: usize) core.errors.Error!void {
    if (index >= count) return error.InvalidAtomIndex;
}

fn countProximity(bonds: []const model.Bond) usize {
    var count: usize = 0;
    for (bonds) |bond| count += @intFromBool(!bond.skip and bond.effective_order == .zero);
    return count;
}

fn flagCrossAtoms(allocator: std.mem.Allocator, atoms: []model.Atom, bonds: []const model.Bond) core.errors.Error!void {
    const degrees = allocator.alloc(u32, atoms.len) catch return error.OutOfMemory;
    defer allocator.free(degrees);
    @memset(degrees, 0);
    for (bonds) |bond| {
        if (bond.skip or bond.effective_order == .zero or atoms[bond.start.index()].hidden or atoms[bond.end.index()].hidden) continue;
        degrees[bond.start.index()] += 1;
        degrees[bond.end.index()] += 1;
    }
    for (atoms) |*atom| atom.cross_layout = !atom.hidden and (atom.atomic_number == .sulfur or atom.atomic_number == .phosphorus);
    for (atoms) |*atom| {
        if (atom.cross_layout or atom.hidden) continue;
        var cross_neighbors: u32 = 0;
        for (bonds) |bond| {
            if (bond.skip or bond.effective_order == .zero or atoms[bond.start.index()].hidden or atoms[bond.end.index()].hidden) continue;
            const neighbor = if (bond.start == atom.id) bond.end else if (bond.end == atom.id) bond.start else continue;
            cross_neighbors += @intFromBool(degrees[neighbor.index()] > 3);
        }
        atom.cross_layout = cross_neighbors > 2;
    }
}

pub const TestAtomStereo = struct {
    value: core.chemistry.AtomStereo = .unspecified,
    looking_from: ?u32 = null,
    atom_a: ?u32 = null,
    atom_b: ?u32 = null,
};

pub const TestAtom = struct {
    atomic_number: core.chemistry.AtomicNumber = .carbon,
    formal_charge: i32 = 0,
    fixed: bool = false,
    constrained: bool = false,
    hidden: bool = false,
    template_coordinates: ?core.math.Vec2 = null,
    coordinates_3d: ?core.math.Vec3 = null,
    stereo: TestAtomStereo = .{},
};

pub const TestBondStereo = struct {
    value: core.chemistry.BondStereo = .unspecified,
    atom_a: ?u32 = null,
    atom_b: ?u32 = null,
};

pub const TestBond = struct {
    start: u32,
    end: u32,
    order: core.chemistry.BondOrder = .single,
    skip: bool = false,
    stereo: TestBondStereo = .{},
    display: core.chemistry.BondDisplay = .none,
    crossing_penalty_multiplier: f32 = 1,
};

pub const TestResidue = struct {
    atom: u32,
    chain: []const u8 = "",
    residue_number: i32 = 0,
    closest_ligand_atom: ?u32 = null,
};

pub const TestInteraction = struct {
    start: u32,
    end: u32,
    other_start_atoms: []const u32 = &.{},
    other_end_atoms: []const u32 = &.{},
    crossing_penalty_multiplier: f32 = 1,
};

pub const TestInput = struct {
    atoms: []const TestAtom,
    bonds: []const TestBond,
    residues: []const TestResidue = &.{},
    residue_interactions: []const TestInteraction = &.{},
    extra_bonds: []const TestBond = &.{},
    options: struct { treat_nonterminal_bonds_to_metal_as_zero_order: bool = true } = .{},
};

const test_atoms = [_]TestAtom{
    .{ .atomic_number = .iron },
    .{ .atomic_number = .carbon },
    .{ .atomic_number = .oxygen },
    .{ .atomic_number = .phosphorus },
    .{ .atomic_number = .carbon, .hidden = true },
};

const test_bonds = [_]TestBond{
    .{ .start = 0, .end = 1 },
    .{ .start = 0, .end = 2 },
    .{ .start = 1, .end = 3 },
    .{ .start = 3, .end = 4 },
};

test "preparation copies input, rewrites nonterminal metal bonds, and preserves structural components" {
    const atoms_before = test_atoms;
    const bonds_before = test_bonds;
    var prepared = try init(std.testing.allocator, TestInput{ .atoms = &test_atoms, .bonds = &test_bonds });
    defer prepared.deinit();

    try std.testing.expectEqualSlices(TestAtom, &atoms_before, &test_atoms);
    try std.testing.expectEqualSlices(TestBond, &bonds_before, &test_bonds);
    try std.testing.expectEqual(@as(u32, 4), prepared.working.active_atom_count);
    try std.testing.expect(prepared.working.atoms[prepared.working.order.input_to_internal[4]].hidden);
    try std.testing.expect(prepared.working.atoms[prepared.working.order.input_to_internal[3]].cross_layout);
    try std.testing.expectEqualSlices(
        u32,
        &.{ 0, 2, 3, 1 },
        prepared.working.order.activeInternalToInput(prepared.working.active_atom_count),
    );
    var active_input_to_internal: [test_atoms.len]u32 = undefined;
    try prepared.working.order.writeActiveInputToInternal(
        prepared.working.active_atom_count,
        &active_input_to_internal,
    );
    try std.testing.expectEqualSlices(
        u32,
        &.{ 0, 3, 1, 2, std.math.maxInt(u32) },
        &active_input_to_internal,
    );

    var effective_by_input: [test_bonds.len]core.chemistry.BondOrder = undefined;
    for (prepared.working.bonds) |bond| effective_by_input[bond.input_index] = bond.effective_order;
    try std.testing.expectEqualSlices(core.chemistry.BondOrder, &.{ .zero, .single, .single, .single }, &effective_by_input);
    try std.testing.expectEqual(@as(usize, 1), prepared.working.proximity_relations.len);
}

test "preparation rejects malformed endpoints and cleans every allocation failure" {
    const invalid_bonds = [_]TestBond{.{ .start = 0, .end = 9 }};
    try std.testing.expectError(error.InvalidAtomIndex, init(
        std.testing.allocator,
        TestInput{ .atoms = &test_atoms, .bonds = &invalid_bonds },
    ));

    const Runner = struct {
        fn run(allocator: std.mem.Allocator, atoms: []const TestAtom, bonds: []const TestBond) !void {
            var prepared = try init(allocator, TestInput{ .atoms = atoms, .bonds = bonds });
            prepared.deinit();
        }
    };
    try core.oom.checkAllocationFailures(std.testing.allocator, Runner.run, .{ &test_atoms, &test_bonds });
}

test "preparation owns residues, extra bonds, interaction ranges, and proximity records" {
    const extra_bonds = [_]TestBond{.{
        .start = 1,
        .end = 2,
        .order = .zero,
        .display = .hashed_reverse,
    }};
    const residues = [_]TestResidue{.{
        .atom = 3,
        .chain = "A",
        .residue_number = 42,
        .closest_ligand_atom = 1,
    }};
    const other_start = [_]u32{ 0, 2 };
    const other_end = [_]u32{1};
    const interactions = [_]TestInteraction{.{
        .start = 3,
        .end = 3,
        .other_start_atoms = &other_start,
        .other_end_atoms = &other_end,
        .crossing_penalty_multiplier = 2,
    }};
    var prepared = try init(std.testing.allocator, TestInput{
        .atoms = &test_atoms,
        .bonds = &test_bonds,
        .extra_bonds = &extra_bonds,
        .residues = &residues,
        .residue_interactions = &interactions,
    });
    defer prepared.deinit();

    try std.testing.expectEqualStrings("A", prepared.working.string_bytes);
    try std.testing.expectEqual(@as(i32, 42), prepared.working.residues[0].residue_number);
    try std.testing.expectEqual(.hashed_reverse, prepared.working.extra_bonds[0].display);
    try std.testing.expectEqual(@as(u32, 2), prepared.working.residue_interactions[0].other_start_atoms.len);
    try std.testing.expectEqual(@as(u32, 1), prepared.working.residue_interactions[0].other_end_atoms.len);
    try std.testing.expectEqual(@as(usize, 3), prepared.working.proximity_relations.len);
    try std.testing.expect(prepared.working.proximity_relations[0] == .bond);
    try std.testing.expect(prepared.working.proximity_relations[1] == .extra_bond);
    try std.testing.expect(prepared.working.proximity_relations[2] == .residue_interaction);
}
