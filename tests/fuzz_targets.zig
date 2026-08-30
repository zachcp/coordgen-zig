//! Fuzz targets for the native surface (cgz-7v2.4.4).
//!
//! Each target runs in two modes without changing shape, which is what makes a
//! promoted seed worth promoting:
//!
//!   `zig build fuzz`   the fuzzer searches for new inputs
//!   `zig build test`   `std.testing.fuzz` replays `options.corpus` and
//!                      nothing else, so every committed seed is a
//!                      deterministic regression test on the ordinary gate
//!
//! Committed seeds live in `tests/fuzz_seeds/` and are embedded below. They are
//! promoted by `tools/run-fuzz`, which proves a seed reproduces before it is
//! written - a seed that reproduces nothing is a regression test that tests
//! nothing, which is cgz-r27's whole subject.
//!
//! ## What the targets assert
//!
//! Not merely "does not crash". A generator that returned zeroes for
//! everything would never crash. Each target asserts the properties the public
//! contract promises, so a silent wrong answer fails too:
//!
//!   * every published observable has the length the input implies
//!   * coordinates are finite
//!   * the two index maps are mutual inverses
//!   * the caller's input is not mutated
//!   * the same input twice produces the same bytes
//!
//! Leak detection is free here: the fuzz runner installs a checking allocator
//! per iteration and exits non-zero if the iteration leaks.

const std = @import("std");
const api = @import("api");
const core = @import("core");

/// Small on purpose. Throughput is the budget: a target that spends its
/// iterations laying out 60-atom molecules explores far less of the decision
/// space than one that runs many small ones, and the interesting branches -
/// ring perception, fragmentation, macrocycle dispatch, component arrangement
/// - are all reachable well below that.
const max_atoms = 12;
const max_bonds = 18;

const elements = [_]core.chemistry.AtomicNumber{
    .carbon, .nitrogen, .oxygen, .sulfur, .fluorine, .chlorine, .phosphorus, .iron, .zinc,
};

const Molecule = struct {
    atom_storage: [max_atoms]api.AtomInput = undefined,
    bond_storage: [max_bonds]api.BondInput = undefined,
    atom_count: usize = 0,
    bond_count: usize = 0,

    fn input(self: *const Molecule) api.Input {
        return .{
            .atoms = self.atom_storage[0..self.atom_count],
            .bonds = self.bond_storage[0..self.bond_count],
        };
    }
};

/// A structurally valid molecule: indices in range, elements ordinary. This
/// target is about the algorithm, not the validator, so it deliberately does
/// not spend iterations on inputs `validate()` rejects at the door -
/// `buildHostile` owns those.
fn buildValid(smith: *std.testing.Smith) Molecule {
    var molecule: Molecule = .{};
    const count = smith.valueRangeAtMost(u8, 1, max_atoms);
    molecule.atom_count = count;
    for (molecule.atom_storage[0..count]) |*atom| {
        atom.* = .{
            .atomic_number = elements[smith.index(elements.len)],
            .formal_charge = smith.valueRangeAtMost(i8, -1, 1),
        };
    }
    var bonds: usize = 0;
    while (bonds < max_bonds and !smith.eosWeightedSimple(6, 1)) {
        const start = smith.index(count);
        const end = smith.index(count);
        if (start == end) continue;
        molecule.bond_storage[bonds] = .{
            .start = @intCast(start),
            .end = @intCast(end),
            // Zero order is a proximity relation rather than an edge, and is
            // drawn rarely so it does not crowd out ordinary connectivity.
            .order = switch (smith.valueRangeAtMost(u8, 0, 9)) {
                0 => .zero,
                1, 2 => .double,
                3 => .triple,
                else => .single,
            },
        };
        bonds += 1;
    }
    molecule.bond_count = bonds;
    return molecule;
}

/// Everything hostile that is REPRESENTABLE: out-of-range endpoints, empty
/// graphs, exotic and virtual elements, extreme charges, contradictory flags.
/// The contract is that these are rejected cleanly or answered well-formed -
/// never that they crash.
///
/// Deliberately no invalid enum values. The Zig API is typed, so an
/// out-of-range `AtomicNumber` or `BondOrder` cannot be passed to it at all;
/// manufacturing one with `@enumFromInt` does not test the library, it just
/// makes the target itself illegal - which is exactly what the first version
/// of this function did, and the fuzzer caught it in 26 runs. Invalid enum
/// values ARE reachable across the C ABI, whose fields are `u32`, and belong
/// to that surface's target. This is the same reason FUZZING.md gives the C
/// ABI its own targets rather than treating it as a thin wrapper.
fn buildHostile(smith: *std.testing.Smith) Molecule {
    var molecule: Molecule = .{};
    const count = smith.valueRangeAtMost(u8, 0, max_atoms);
    molecule.atom_count = count;
    for (molecule.atom_storage[0..count]) |*atom| {
        atom.* = .{
            // Every element the enum names, 0 (virtual) through 118, not just
            // the organic handful the valid target draws from.
            .atomic_number = @fromBackingInt(@intCast(smith.valueRangeAtMost(u8, 0, 118))),
            .formal_charge = smith.value(i16),
            .hidden = smith.boolWeighted(9, 1),
            .fixed = smith.boolWeighted(9, 1),
            .constrained = smith.boolWeighted(9, 1),
        };
    }
    var bonds: usize = 0;
    while (bonds < max_bonds and !smith.eosWeightedSimple(6, 1)) {
        molecule.bond_storage[bonds] = .{
            // Unconstrained by `count`, so out-of-range endpoints - the main
            // hostility the type system still permits - are reached often.
            .start = smith.value(u8),
            .end = smith.value(u8),
            .order = switch (smith.valueRangeAtMost(u8, 0, 3)) {
                0 => .zero,
                1 => .single,
                2 => .double,
                else => .triple,
            },
            .skip = smith.boolWeighted(9, 1),
        };
        bonds += 1;
    }
    molecule.bond_count = bonds;
    return molecule;
}

/// Every property the public result promises. Anything a generator returning
/// zeroes would satisfy is not worth asserting here.
fn checkResult(molecule: *const Molecule, result: api.Result) !void {
    const input = molecule.input();
    try std.testing.expectEqual(input.atoms.len, result.coordinates.len);
    try std.testing.expectEqual(input.atoms.len, result.input_to_internal.len);
    try std.testing.expectEqual(input.atoms.len, result.internal_to_input.len);
    try std.testing.expectEqual(input.atoms.len, result.atom_stereo.len);
    try std.testing.expectEqual(input.bonds.len, result.effective_bond_orders.len);
    try std.testing.expectEqual(input.bonds.len, result.bond_displays.len);

    for (result.coordinates) |coordinate| {
        try std.testing.expect(coordinate.isFinite());
    }
    // Mutual inverses. A map that merely stays in range would pass a bounds
    // check while scrambling the caller's atom identities.
    for (result.input_to_internal, 0..) |internal, caller_index| {
        try std.testing.expect(internal < result.internal_to_input.len);
        try std.testing.expectEqual(caller_index, result.internal_to_input[internal]);
    }
}

fn generateOnce(molecule: *const Molecule) !?api.Result {
    return api.generate(std.testing.allocator, molecule.input()) catch |err| switch (err) {
        // Enumerated rather than caught wholesale, and deliberately with no
        // `else`: the error set is closed, so adding a variant to the public
        // API stops this compiling until someone decides whether the fuzzer
        // should treat it as a refusal or as a finding. A compile error is a
        // better guarantee than a runtime check, and this is the only place
        // that gets one for free.
        error.EmptyGraph,
        error.TooManyItems,
        error.InvalidAtomicNumber,
        error.InvalidBondOrder,
        error.InvalidAtomIndex,
        error.InvalidStereo,
        error.InvalidCoordinate,
        error.InvalidOption,
        error.InvalidMapping,
        error.Unsupported,
        => return null,
        error.OutOfMemory => return null,
    };
}

fn generationInvariants(_: void, smith: *std.testing.Smith) anyerror!void {
    var molecule = buildValid(smith);
    const before = molecule;

    var first = (try generateOnce(&molecule)) orelse return;
    defer first.deinit();
    try checkResult(&molecule, first);

    // The caller's input is borrowed, never written through.
    try std.testing.expectEqualSlices(
        api.AtomInput,
        before.atom_storage[0..before.atom_count],
        molecule.atom_storage[0..molecule.atom_count],
    );
    try std.testing.expectEqualSlices(
        api.BondInput,
        before.bond_storage[0..before.bond_count],
        molecule.bond_storage[0..molecule.bond_count],
    );

    // Determinism, on every input the fuzzer reaches rather than on the ten
    // families tests/native_determinism.zig can enumerate by hand.
    var second = (try generateOnce(&molecule)) orelse return;
    defer second.deinit();
    try std.testing.expectEqual(first.clean_pose, second.clean_pose);
    try std.testing.expectEqualSlices(api.Vec2, first.coordinates, second.coordinates);
    try std.testing.expectEqualSlices(u32, first.input_to_internal, second.input_to_internal);
    try std.testing.expectEqualSlices(api.BondOrder, first.effective_bond_orders, second.effective_bond_orders);
}

fn hostileInputIsRejected(_: void, smith: *std.testing.Smith) anyerror!void {
    var molecule = buildHostile(smith);
    var result = (try generateOnce(&molecule)) orelse return;
    defer result.deinit();
    // A hostile input that IS accepted has to be as well-formed as any other
    // result. Accepting it is allowed; accepting it and returning nonsense is
    // not.
    try checkResult(&molecule, result);
}

const domain_max_atoms = 16;
const domain_max_bonds = 16;

const DomainInput = struct {
    atoms: [domain_max_atoms]api.AtomInput = undefined,
    bonds: [domain_max_bonds]api.BondInput = undefined,
    residues: [4]api.ResidueInput = undefined,
    interactions: [2]api.ResidueInteractionInput = undefined,
    atom_count: usize = 0,
    bond_count: usize = 0,
    residue_count: usize = 0,
    interaction_count: usize = 0,
    options: api.Options = .{},

    fn input(self: *const DomainInput) api.Input {
        return .{
            .atoms = self.atoms[0..self.atom_count],
            .bonds = self.bonds[0..self.bond_count],
            .residues = self.residues[0..self.residue_count],
            .residue_interactions = self.interactions[0..self.interaction_count],
            .options = self.options,
        };
    }
};

/// Reach domains the small graph targets cannot: residue/protein placement,
/// the 13-member macrocycle dispatch, and every public option. Unsupported
/// runtime template directories are intentional inputs too: until a native
/// loader exists, explicit rejection is the public contract to preserve.
fn buildDomainInput(smith: *std.testing.Smith) DomainInput {
    var built: DomainInput = .{};
    switch (smith.valueRangeAtMost(u8, 0, 2)) {
        0 => {
            built.atom_count = 4;
            for (built.atoms[0..built.atom_count]) |*atom| atom.* = .{};
            built.bonds[0] = .{ .start = 0, .end = 1 };
            built.bonds[1] = .{ .start = 1, .end = 2 };
            built.bond_count = 2;
            built.residues[0] = .{
                .atom = 3,
                .chain = "A",
                .residue_number = smith.valueRangeAtMost(i8, -4, 4),
                .closest_ligand_atom = 1,
            };
            built.residue_count = 1;
            built.interactions[0] = .{
                .start = 3,
                .end = 1,
                .crossing_penalty_multiplier = smith.valueRangeAtMost(u8, 0, 4),
            };
            built.interaction_count = 1;
        },
        1 => {
            built.atom_count = 4;
            for (built.atoms[0..built.atom_count]) |*atom| atom.* = .{};
            built.residues[0] = .{ .atom = 0, .chain = "A", .residue_number = 1 };
            built.residues[1] = .{ .atom = 1, .chain = "A", .residue_number = 2 };
            built.residues[2] = .{ .atom = 2, .chain = "B", .residue_number = 1 };
            built.residues[3] = .{ .atom = 3, .chain = "B", .residue_number = 2 };
            built.residue_count = 4;
            built.interactions[0] = .{ .start = 0, .end = 2 };
            built.interactions[1] = .{ .start = 1, .end = 3 };
            built.interaction_count = 2;
        },
        else => {
            built.atom_count = 13;
            built.bond_count = 13;
            for (built.atoms[0..built.atom_count], built.bonds[0..built.bond_count], 0..) |*atom, *bond, index| {
                atom.* = .{};
                bond.* = .{
                    .start = @intCast(index),
                    .end = @intCast((index + 1) % built.atom_count),
                };
            }
        },
    }

    built.options = .{
        .precision = switch (smith.valueRangeAtMost(u8, 0, 2)) {
            0 => api.Precision.quick,
            1 => api.Precision.standard,
            else => api.Precision.best,
        },
        .score_residue_interactions = smith.boolWeighted(1, 1),
        .treat_nonterminal_bonds_to_metal_as_zero_order = smith.boolWeighted(1, 1),
        .even_angles = smith.boolWeighted(15, 1),
        .skip_minimization = smith.boolWeighted(1, 1),
        .force_open_macrocycles = smith.boolWeighted(1, 1),
        .constrain_all_atoms = smith.boolWeighted(15, 1),
        .debug_coordinates = smith.boolWeighted(15, 1),
        .load_templates = smith.boolWeighted(1, 1),
        .template_directory = if (smith.boolWeighted(31, 1)) "fixtures" else null,
    };
    return built;
}

fn domainResult(input: api.Input) !?api.Result {
    return api.generate(std.testing.allocator, input) catch |err| switch (err) {
        error.EmptyGraph,
        error.TooManyItems,
        error.InvalidAtomicNumber,
        error.InvalidBondOrder,
        error.InvalidAtomIndex,
        error.InvalidStereo,
        error.InvalidCoordinate,
        error.InvalidOption,
        error.InvalidMapping,
        error.Unsupported,
        error.OutOfMemory,
        => return null,
    };
}

fn residueTemplateOptionsAndMacrocycle(_: void, smith: *std.testing.Smith) anyerror!void {
    var built = buildDomainInput(smith);
    const input = built.input();
    var first = (try domainResult(input)) orelse return;
    defer first.deinit();
    try std.testing.expectEqual(input.atoms.len, first.coordinates.len);
    for (first.coordinates) |coordinate| try std.testing.expect(coordinate.isFinite());

    var second = (try domainResult(input)) orelse return error.NondeterministicDomainAcceptance;
    defer second.deinit();
    try std.testing.expectEqual(first.clean_pose, second.clean_pose);
    try std.testing.expectEqualSlices(api.Vec2, first.coordinates, second.coordinates);
    try std.testing.expectEqualSlices(u32, first.input_to_internal, second.input_to_internal);
}

/// Seeds promoted by tools/run-fuzz. Each was proven to reproduce its failure
/// before being written here, and each runs on ordinary `zig build test`.
const generation_seeds = @import("fuzz_seeds/generation_invariants_hold_on_structurally_valid_molecules.zig").seeds;
const hostile_seeds = @import("fuzz_seeds/hostile_input_is_rejected_or_answered_well_formed.zig").seeds;
const domain_seeds = @import("fuzz_seeds/residue_template_options_and_macrocycle_domains.zig").seeds;

test "fuzz: generation invariants hold on structurally valid molecules" {
    try std.testing.fuzz({}, generationInvariants, .{ .corpus = generation_seeds });
}

test "fuzz: hostile input is rejected or answered well-formed" {
    try std.testing.fuzz({}, hostileInputIsRejected, .{ .corpus = hostile_seeds });
}

test "fuzz: residue template options and macrocycle domains" {
    try std.testing.fuzz({}, residueTemplateOptionsAndMacrocycle, .{ .corpus = domain_seeds });
}
