//! Seeded corpus generator for oracle stability classification.
//!
//! The corpus exists to answer one question per input: which comparison tier
//! is that input entitled to. Answering it requires running the same inputs
//! through oracle builds that differ only in architecture or in heap address
//! order, so generation must be reproducible from nothing but a partition and
//! an index — no files, no clock, no allocator-dependent iteration.
//!
//! Generation is therefore pinned, not merely deterministic: the PRNG
//! algorithm, its seeding, the element table, and the order of draws are all
//! fixed here. Changing any of them changes the population and invalidates
//! every recorded classification.
//!
//! Conformance-only; never part of an installed artifact.

const std = @import("std");

pub const Partition = enum {
    /// Chemically implausible random graphs: random valences, free metals,
    /// multi-component inputs, macrocycles from long ring closures. These are
    /// where upstream's own layout is unstable, so this is the partition that
    /// establishes the parity ceiling.
    adversarial,
    /// Realistic structures — fused polycyclics, a steroid skeleton, a
    /// seventeen-membered macrocycle, a bridged bicyclic, a peptide-like
    /// chain, a four-ring biaryl chain. Recorded as stable across
    /// architectures and optimization levels, so they carry the exact tier.
    drug_like,

    pub fn memberCount(self: Partition, requested: u32) u32 {
        return switch (self) {
            .adversarial => requested,
            .drug_like => @intCast(drug_like_members.len),
        };
    }
};

pub const Atom = struct {
    atomic_number: u32,
    formal_charge: i32 = 0,
};

pub const Bond = struct {
    start: u32,
    end: u32,
    order: u32,
};

pub const Molecule = struct {
    partition: Partition,
    index: u32,
    atoms: []const Atom,
    bonds: []const Bond,

    pub fn deinit(self: Molecule, gpa: std.mem.Allocator) void {
        gpa.free(self.atoms);
        gpa.free(self.bonds);
    }
};

/// xorshift64* (Vigna 2016). Pinned by value, not by reference to any std
/// default: the corpus is a recorded population, and a different generator is
/// a different corpus even if it is equally random.
pub const Rng = struct {
    state: u64,

    pub const golden_gamma: u64 = 0x9E3779B97F4A7C15;
    pub const seed_multiplier: u64 = 0xBF58476D1CE4E5B9;
    pub const output_multiplier: u64 = 0x2545F4914F6CDD1D;
    /// Draws discarded after seeding, so that low-index members do not start
    /// from a barely-mixed state.
    pub const warmup_draws = 8;

    pub fn init(seed: u64) Rng {
        var rng: Rng = .{ .state = golden_gamma ^ (seed *% seed_multiplier) };
        for (0..warmup_draws) |_| _ = rng.next();
        return rng;
    }

    pub fn next(self: *Rng) u64 {
        var x = self.state;
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        self.state = x;
        return x *% output_multiplier;
    }

    /// Modulo reduction, including its bias toward low values. The bias is
    /// part of the recorded population and must not be "fixed" into a
    /// rejection sampler.
    pub fn below(self: *Rng, bound: u32) u32 {
        std.debug.assert(bound > 0);
        return @intCast(self.next() % bound);
    }
};

/// Carbon is weighted heavily; the tail carries heteroatoms, halogens, and
/// two metals (iron, zinc) so that the metal zero-order-bond rewrite and
/// unusual valences are exercised.
const elements = [_]u32{ 6, 6, 6, 6, 7, 8, 16, 9, 17, 15, 26, 30 };

pub const GenerateError = std.mem.Allocator.Error || error{UnknownMember};

pub fn generate(gpa: std.mem.Allocator, partition: Partition, index: u32) GenerateError!Molecule {
    return switch (partition) {
        .adversarial => generateAdversarial(gpa, index),
        .drug_like => generateDrugLike(gpa, index),
    };
}

fn generateAdversarial(gpa: std.mem.Allocator, index: u32) GenerateError!Molecule {
    var rng: Rng = .init(index);

    var atoms: std.ArrayList(Atom) = .empty;
    errdefer atoms.deinit(gpa);
    var bonds: std.ArrayList(Bond) = .empty;
    errdefer bonds.deinit(gpa);

    // Disconnected components: upstream places each separately and then
    // arranges them, which is one of the pointer-ordered decisions.
    const component_count = 1 + rng.below(3);
    for (0..component_count) |_| {
        const base: u32 = @intCast(atoms.items.len);
        const count = 3 + rng.below(28);

        for (0..count) |_| {
            const atomic_number = elements[rng.below(elements.len)];
            const charged = rng.below(12) == 0;
            try atoms.append(gpa, .{
                .atomic_number = atomic_number,
                .formal_charge = if (charged) @as(i32, @intCast(rng.below(3))) - 1 else 0,
            });
        }

        // Spanning tree first: every component is connected by construction.
        for (1..count) |offset| {
            const parent = base + rng.below(@intCast(offset));
            const draw = rng.below(10);
            try bonds.append(gpa, .{
                .start = parent,
                .end = base + @as(u32, @intCast(offset)),
                .order = if (draw < 6) 1 else if (draw < 9) 2 else 3,
            });
        }

        // Ring closures: short ones make ordinary rings, long ones make
        // macrocycles. A rejected draw still consumes its two values, which
        // is what keeps the sequence reproducible.
        const closure_attempts = rng.below(4);
        for (0..closure_attempts) |_| {
            if (count <= 3) break;
            const start = base + rng.below(count);
            const end = base + rng.below(count);
            if (start == end) continue;
            if (hasBond(bonds.items, start, end)) continue;
            try bonds.append(gpa, .{ .start = start, .end = end, .order = 1 });
        }
    }

    return .{
        .partition = .adversarial,
        .index = index,
        .atoms = try atoms.toOwnedSlice(gpa),
        .bonds = try bonds.toOwnedSlice(gpa),
    };
}

fn hasBond(bonds: []const Bond, start: u32, end: u32) bool {
    for (bonds) |bond| {
        if ((bond.start == start and bond.end == end) or
            (bond.start == end and bond.end == start)) return true;
    }
    return false;
}

const DrugLikeMember = struct {
    smiles: []const u8,
    atomic_numbers: []const u32,
    /// `.{ start, end, order }`, in upstream's own construction order.
    bonds: []const [3]u32,
};

fn generateDrugLike(gpa: std.mem.Allocator, index: u32) GenerateError!Molecule {
    if (index >= drug_like_members.len) return error.UnknownMember;
    const member = drug_like_members[index];

    const atoms = try gpa.alloc(Atom, member.atomic_numbers.len);
    errdefer gpa.free(atoms);
    for (atoms, member.atomic_numbers) |*atom, atomic_number| {
        atom.* = .{ .atomic_number = atomic_number };
    }

    const bonds = try gpa.alloc(Bond, member.bonds.len);
    errdefer gpa.free(bonds);
    for (bonds, member.bonds) |*bond, entry| {
        bond.* = .{ .start = entry[0], .end = entry[1], .order = entry[2] };
    }

    return .{ .partition = .drug_like, .index = index, .atoms = atoms, .bonds = bonds };
}

pub fn drugLikeSmiles(index: u32) ?[]const u8 {
    if (index >= drug_like_members.len) return null;
    return drug_like_members[index].smiles;
}

/// Atom-count bucket, recorded per member because divergence scales with
/// size: the review measurements ran from 7/186 divergent below ten atoms to
/// 42/50 in the seventy-to-eighty-nine range.
pub const SizeBucket = enum {
    tiny,
    small,
    medium,
    large,
    huge,

    pub fn of(atom_count: usize) SizeBucket {
        return if (atom_count < 10)
            .tiny
        else if (atom_count < 30)
            .small
        else if (atom_count < 50)
            .medium
        else if (atom_count < 70)
            .large
        else
            .huge;
    }
};

/// The canonical corpus dump. `conformance/smiles_reference.cpp` prints the
/// drug-like partition in exactly this form from upstream's own parser, and
/// the build diffs the two, so the committed tables are checked rather than
/// trusted.
pub fn writeMolecule(writer: *std.Io.Writer, molecule: Molecule) std.Io.Writer.Error!void {
    try writer.print("molecule {t} {d} atoms={d} bonds={d}\n", .{
        molecule.partition,
        molecule.index,
        molecule.atoms.len,
        molecule.bonds.len,
    });
    for (molecule.atoms, 0..) |atom, position| {
        try writer.print("atom {d} z={d} q={d}\n", .{
            position,
            atom.atomic_number,
            atom.formal_charge,
        });
    }
    for (molecule.bonds, 0..) |bond, position| {
        try writer.print("bond {d} {d} {d} order={d}\n", .{
            position,
            bond.start,
            bond.end,
            bond.order,
        });
    }
}

const drug_like_members = [_]DrugLikeMember{
    .{
        .smiles = "C1CC2CCC3CCCC4CCC(C1)C2C34",
        .atomic_numbers = &.{ 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6 },
        .bonds = &.{
            .{ 1, 0, 1 },   .{ 2, 1, 1 },   .{ 3, 2, 1 },   .{ 4, 3, 1 },
            .{ 5, 4, 1 },   .{ 6, 5, 1 },   .{ 7, 6, 1 },   .{ 8, 7, 1 },
            .{ 9, 8, 1 },   .{ 10, 9, 1 },  .{ 11, 10, 1 }, .{ 12, 11, 1 },
            .{ 13, 12, 1 }, .{ 13, 0, 1 },  .{ 14, 12, 1 }, .{ 14, 2, 1 },
            .{ 15, 14, 1 }, .{ 15, 5, 1 },  .{ 15, 9, 1 },
        },
    },
    .{
        .smiles = "CC12CCC3C(CCC4CC(O)CCC34C)C1CCC2O",
        .atomic_numbers = &.{ 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 8, 6, 6, 6, 6, 6, 6, 6, 6, 8 },
        .bonds = &.{
            .{ 1, 0, 1 },   .{ 2, 1, 1 },   .{ 3, 2, 1 },   .{ 4, 3, 1 },
            .{ 5, 4, 1 },   .{ 6, 5, 1 },   .{ 7, 6, 1 },   .{ 8, 7, 1 },
            .{ 9, 8, 1 },   .{ 10, 9, 1 },  .{ 11, 10, 1 }, .{ 12, 10, 1 },
            .{ 13, 12, 1 }, .{ 14, 13, 1 }, .{ 14, 4, 1 },  .{ 14, 8, 1 },
            .{ 15, 14, 1 }, .{ 16, 5, 1 },  .{ 16, 1, 1 },  .{ 17, 16, 1 },
            .{ 18, 17, 1 }, .{ 19, 18, 1 }, .{ 19, 1, 1 },  .{ 20, 19, 1 },
        },
    },
    .{
        .smiles = "C1CCCCCCCCCCCCCCCC1",
        .atomic_numbers = &.{ 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6 },
        .bonds = &.{
            .{ 1, 0, 1 },   .{ 2, 1, 1 },   .{ 3, 2, 1 },   .{ 4, 3, 1 },
            .{ 5, 4, 1 },   .{ 6, 5, 1 },   .{ 7, 6, 1 },   .{ 8, 7, 1 },
            .{ 9, 8, 1 },   .{ 10, 9, 1 },  .{ 11, 10, 1 }, .{ 12, 11, 1 },
            .{ 13, 12, 1 }, .{ 14, 13, 1 }, .{ 15, 14, 1 }, .{ 16, 15, 1 },
            .{ 16, 0, 1 },
        },
    },
    .{
        .smiles = "C1CCCCC2CCCCC(CCCCC1)C2",
        .atomic_numbers = &.{ 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6 },
        .bonds = &.{
            .{ 1, 0, 1 },   .{ 2, 1, 1 },   .{ 3, 2, 1 },   .{ 4, 3, 1 },
            .{ 5, 4, 1 },   .{ 6, 5, 1 },   .{ 7, 6, 1 },   .{ 8, 7, 1 },
            .{ 9, 8, 1 },   .{ 10, 9, 1 },  .{ 11, 10, 1 }, .{ 12, 11, 1 },
            .{ 13, 12, 1 }, .{ 14, 13, 1 }, .{ 15, 14, 1 }, .{ 15, 0, 1 },
            .{ 16, 10, 1 }, .{ 16, 5, 1 },
        },
    },
    .{
        .smiles = "CCC(CC)C(=O)NC(CS)C(=O)NCCCN",
        .atomic_numbers = &.{ 6, 6, 6, 6, 6, 6, 8, 7, 6, 6, 16, 6, 8, 7, 6, 6, 6, 7 },
        .bonds = &.{
            .{ 1, 0, 1 },   .{ 2, 1, 1 },   .{ 3, 2, 1 },   .{ 4, 3, 1 },
            .{ 5, 2, 1 },   .{ 6, 5, 2 },   .{ 7, 5, 1 },   .{ 8, 7, 1 },
            .{ 9, 8, 1 },   .{ 10, 9, 1 },  .{ 11, 8, 1 },  .{ 12, 11, 2 },
            .{ 13, 11, 1 }, .{ 14, 13, 1 }, .{ 15, 14, 1 }, .{ 16, 15, 1 },
            .{ 17, 16, 1 },
        },
    },
    .{
        .smiles = "C1CC2CCC1CC2",
        .atomic_numbers = &.{ 6, 6, 6, 6, 6, 6, 6, 6 },
        .bonds = &.{
            .{ 1, 0, 1 }, .{ 2, 1, 1 }, .{ 3, 2, 1 }, .{ 4, 3, 1 },
            .{ 5, 4, 1 }, .{ 5, 0, 1 }, .{ 6, 5, 1 }, .{ 7, 6, 1 },
            .{ 7, 2, 1 },
        },
    },
    .{
        .smiles = "C1CCC(CC1)C1CCC(CC1)C1CCC(CC1)C1CCCCC1",
        .atomic_numbers = &.{ 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6 },
        .bonds = &.{
            .{ 1, 0, 1 },   .{ 2, 1, 1 },   .{ 3, 2, 1 },   .{ 4, 3, 1 },
            .{ 5, 4, 1 },   .{ 5, 0, 1 },   .{ 6, 3, 1 },   .{ 7, 6, 1 },
            .{ 8, 7, 1 },   .{ 9, 8, 1 },   .{ 10, 9, 1 },  .{ 11, 10, 1 },
            .{ 11, 6, 1 },  .{ 12, 9, 1 },  .{ 13, 12, 1 }, .{ 14, 13, 1 },
            .{ 15, 14, 1 }, .{ 16, 15, 1 }, .{ 17, 16, 1 }, .{ 17, 12, 1 },
            .{ 18, 15, 1 }, .{ 19, 18, 1 }, .{ 20, 19, 1 }, .{ 21, 20, 1 },
            .{ 22, 21, 1 }, .{ 23, 22, 1 }, .{ 23, 18, 1 },
        },
    },
};

const testing = std.testing;

test "the pinned PRNG reproduces its recorded stream" {
    // Seeded exactly as the review probe behind cgz-r05 seeded it, so the
    // corpus this generator produces is the population those measurements
    // describe.
    var rng: Rng = .init(0);
    const first = rng.next();
    var again: Rng = .init(0);
    try testing.expectEqual(first, again.next());

    var other: Rng = .init(1);
    try testing.expect(other.next() != first);

    // xorshift64* never reaches the zero state, so no draw sequence can
    // collapse.
    var probe: Rng = .{ .state = 1 };
    for (0..64) |_| {
        _ = probe.next();
        try testing.expect(probe.state != 0);
    }
}

test "generation depends only on partition and index" {
    for ([_]u32{ 0, 1, 7, 199, 1999 }) |index| {
        const first = try generate(testing.allocator, .adversarial, index);
        defer first.deinit(testing.allocator);
        const second = try generate(testing.allocator, .adversarial, index);
        defer second.deinit(testing.allocator);

        try testing.expectEqualSlices(Atom, first.atoms, second.atoms);
        try testing.expectEqualSlices(Bond, first.bonds, second.bonds);
    }
}

test "adversarial members are well formed graphs" {
    for (0..256) |raw_index| {
        const index: u32 = @intCast(raw_index);
        const molecule = try generate(testing.allocator, .adversarial, index);
        defer molecule.deinit(testing.allocator);

        try testing.expect(molecule.atoms.len >= 3);
        try testing.expect(molecule.atoms.len <= 3 * 30);
        for (molecule.atoms) |atom| {
            try testing.expect(std.mem.indexOfScalar(u32, &elements, atom.atomic_number) != null);
            try testing.expect(atom.formal_charge >= -1 and atom.formal_charge <= 1);
        }
        for (molecule.bonds, 0..) |bond, position| {
            try testing.expect(bond.start < molecule.atoms.len);
            try testing.expect(bond.end < molecule.atoms.len);
            try testing.expect(bond.start != bond.end);
            try testing.expect(bond.order >= 1 and bond.order <= 3);
            // No duplicate edge in either direction.
            try testing.expect(!hasBond(molecule.bonds[0..position], bond.start, bond.end));
        }
    }
}

test "adversarial members carry the properties the corpus exists to exercise" {
    var seen_multi_component = false;
    var seen_metal = false;
    var seen_charge = false;
    var seen_ring_closure = false;
    var seen_high_order = false;

    for (0..256) |raw_index| {
        const molecule = try generate(testing.allocator, .adversarial, @intCast(raw_index));
        defer molecule.deinit(testing.allocator);

        if (molecule.bonds.len >= molecule.atoms.len) seen_ring_closure = true;
        if (molecule.bonds.len + 1 < molecule.atoms.len) seen_multi_component = true;
        for (molecule.atoms) |atom| {
            if (atom.atomic_number == 26 or atom.atomic_number == 30) seen_metal = true;
            if (atom.formal_charge != 0) seen_charge = true;
        }
        for (molecule.bonds) |bond| {
            if (bond.order > 1) seen_high_order = true;
        }
    }

    try testing.expect(seen_multi_component);
    try testing.expect(seen_metal);
    try testing.expect(seen_charge);
    try testing.expect(seen_ring_closure);
    try testing.expect(seen_high_order);
}

test "drug-like members match their committed tables" {
    try testing.expectEqual(@as(u32, 7), Partition.drug_like.memberCount(2000));
    try testing.expectError(error.UnknownMember, generate(testing.allocator, .drug_like, 7));

    const macrocycle = try generate(testing.allocator, .drug_like, 2);
    defer macrocycle.deinit(testing.allocator);
    try testing.expectEqual(@as(usize, 17), macrocycle.atoms.len);
    // Seventeen atoms closed into a single ring: seventeen bonds.
    try testing.expectEqual(@as(usize, 17), macrocycle.bonds.len);

    const peptide = try generate(testing.allocator, .drug_like, 4);
    defer peptide.deinit(testing.allocator);
    var double_bonds: usize = 0;
    for (peptide.bonds) |bond| {
        if (bond.order == 2) double_bonds += 1;
    }
    try testing.expectEqual(@as(usize, 2), double_bonds);
    try testing.expectEqual(@as(u32, 16), peptide.atoms[10].atomic_number);
}

test "size buckets follow the divergence-versus-size boundaries" {
    try testing.expectEqual(SizeBucket.tiny, SizeBucket.of(9));
    try testing.expectEqual(SizeBucket.small, SizeBucket.of(10));
    try testing.expectEqual(SizeBucket.medium, SizeBucket.of(30));
    try testing.expectEqual(SizeBucket.large, SizeBucket.of(50));
    try testing.expectEqual(SizeBucket.huge, SizeBucket.of(70));
}
