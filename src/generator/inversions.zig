//! Chiral inversion constraints for the minimization interaction set
//! (cgz-r34).
//!
//! Mirrors `CoordgenMinimizer::addChiralInversionConstraintsOfMolecule`, which
//! is narrower than its name suggests: it walks macrocycle rings only, and for
//! each ring atom whose preceding ring bond carries stereochemistry it emits
//! one E/Z constraint over the four consecutive atoms around that bond. Small
//! rings are excluded because their geometry is already fixed by the ring
//! itself.
//!
//! Lives in the generator layer for the same reason as bend construction: the
//! approved module edges give `optimize` neither topology nor the graph.

const std = @import("std");
const core = @import("core");
const model = @import("model");
const topology = @import("topology");

/// Build the constraints in ring order, then atom order within each ring,
/// which is upstream's emission order.
pub fn buildChiralInversionConstraints(
    allocator: std.mem.Allocator,
    atoms: []const model.Atom,
    bonds: []const model.Bond,
    graph: topology.Graph,
    membership: topology.RingMembership,
) core.errors.Error![]core.interaction.EzConstraint {
    var constraints: std.ArrayList(core.interaction.EzConstraint) = .empty;
    defer constraints.deinit(allocator);

    for (membership.rings) |ring| {
        const ring_atoms = membership.atoms(ring.id);
        if (ring_atoms.len < topology.rings.macrocycle_size) continue;
        const size = ring_atoms.len;
        for (0..size) |index| {
            const previous = ring_atoms[(index + size - 1) % size];
            const before_previous = ring_atoms[(index + size - 2) % size];
            const next = ring_atoms[(index + 1) % size];
            const current = ring_atoms[index];

            const bond_id = bondBetween(graph, previous, current) orelse continue;
            if (bond_id.index() >= bonds.len) return error.InvalidMapping;
            const bond = bonds[bond_id.index()];
            if (!topology.stereo.isStereoBond(membership, bond)) continue;

            const cis = try markedAsCis(
                allocator,
                atoms,
                bonds,
                graph,
                membership,
                bond,
                before_previous,
                next,
            );
            constraints.append(allocator, .{
                .side_a = before_previous,
                .double_a = previous,
                .double_b = current,
                .side_b = next,
                .is_z = cis,
            }) catch return error.OutOfMemory;
        }
    }
    return constraints.toOwnedSlice(allocator) catch return error.OutOfMemory;
}

/// `sketcherMinimizerBond::markedAsCis`: start from the bond's absolute Z
/// state, then invert once for each reference atom that is not a CIP first
/// neighbour of either endpoint. The two inversions cancel when both
/// references are off the CIP pair, which is why this is not simply "is Z".
fn markedAsCis(
    allocator: std.mem.Allocator,
    atoms: []const model.Atom,
    bonds: []const model.Bond,
    graph: topology.Graph,
    membership: topology.RingMembership,
    bond: model.Bond,
    first: core.ids.AtomId,
    second: core.ids.AtomId,
) core.errors.Error!bool {
    const absolute = try topology.stereo.absoluteBondStereo(
        allocator,
        atoms,
        bonds,
        graph,
        membership,
        bond.id,
    );
    var cis = absolute == .z;
    const start_neighbor = try topology.stereo.firstNeighbor(
        allocator,
        atoms,
        bonds,
        graph,
        bond.start,
        bond.end,
    );
    const end_neighbor = try topology.stereo.firstNeighbor(
        allocator,
        atoms,
        bonds,
        graph,
        bond.end,
        bond.start,
    );
    if (!isCipFirst(first, start_neighbor, end_neighbor)) cis = !cis;
    if (!isCipFirst(second, start_neighbor, end_neighbor)) cis = !cis;
    return cis;
}

fn isCipFirst(
    atom: core.ids.AtomId,
    start_neighbor: ?core.ids.AtomId,
    end_neighbor: ?core.ids.AtomId,
) bool {
    if (start_neighbor) |candidate| {
        if (candidate == atom) return true;
    }
    if (end_neighbor) |candidate| {
        if (candidate == atom) return true;
    }
    return false;
}

/// `CoordgenMinimizer::getChetoCs`' membership test for one atom: a carbon
/// carrying a double bond to an oxygen. Shared with the global orientation
/// stage, which classifies peptides from the same three predicates without the
/// two-of-each-class guard this file's constraint builder applies (cgz-r31.1).
pub fn isChetoCarbon(
    atom: core.ids.AtomId,
    atoms: []const model.Atom,
    bonds: []const model.Bond,
    graph: topology.Graph,
) core.errors.Error!bool {
    if (atom.index() >= atoms.len) return error.InvalidMapping;
    if (atoms[atom.index()].atomic_number != .carbon) return false;
    for (graph.neighbors(atom)) |neighbor| {
        if (neighbor.index() >= atoms.len) return error.InvalidMapping;
        if (atoms[neighbor.index()].atomic_number != .oxygen) continue;
        const bond_id = bondBetween(graph, atom, neighbor) orelse continue;
        if (bond_id.index() >= bonds.len) return error.InvalidMapping;
        if (bonds[bond_id.index()].effective_order != .double) continue;
        return true;
    }
    return false;
}

fn bondBetween(
    graph: topology.Graph,
    first: core.ids.AtomId,
    second: core.ids.AtomId,
) ?core.ids.BondId {
    for (graph.incidentBonds(first)) |candidate| {
        for (graph.incidentBonds(second)) |other| {
            if (candidate == other) return candidate;
        }
    }
    return null;
}

const testing = std.testing;

test "only macrocycle stereo bonds produce inversion constraints" {
    // A ten-membered carbocycle with one stereo double bond, plus a small
    // ring case that must produce nothing.
    const ring_size = 10;
    var atoms: [ring_size]model.Atom = undefined;
    var bonds: [ring_size]model.Bond = undefined;
    for (&atoms, &bonds, 0..) |*atom, *bond, index| {
        const angle = @as(f32, @floatFromInt(index)) * std.math.tau / ring_size;
        atom.* = .{
            .id = core.ids.AtomId.fromIndex(@intCast(index)),
            .input_index = @intCast(index),
            .atomic_number = .carbon,
            .coordinates = .{ .x = 100 * @cos(angle), .y = 100 * @sin(angle) },
        };
        bond.* = .{
            .id = core.ids.BondId.fromIndex(@intCast(index)),
            .input_index = @intCast(index),
            .start = core.ids.AtomId.fromIndex(@intCast(index)),
            .end = core.ids.AtomId.fromIndex(@intCast((index + 1) % ring_size)),
            .input_order = .single,
            .effective_order = .single,
        };
    }
    // Make bond 0 a stereo double bond with explicit references.
    bonds[0].input_order = .double;
    bonds[0].effective_order = .double;
    bonds[0].stereo = .cis;
    bonds[0].stereo_atom_a = core.ids.AtomId.fromIndex(ring_size - 1);
    bonds[0].stereo_atom_b = core.ids.AtomId.fromIndex(2);

    var graph = try topology.Graph.init(testing.allocator, &atoms, &bonds);
    defer graph.deinit();
    var membership = try topology.RingMembership.init(testing.allocator, graph, &bonds);
    defer membership.deinit();

    const constraints = try buildChiralInversionConstraints(
        testing.allocator,
        &atoms,
        &bonds,
        graph,
        membership,
    );
    defer testing.allocator.free(constraints);

    // Exactly the one stereo bond contributes, and its constraint spans the
    // four consecutive ring atoms around it.
    try testing.expectEqual(@as(usize, 1), constraints.len);
    try testing.expectEqual(core.ids.AtomId.fromIndex(0), constraints[0].double_a);
    try testing.expectEqual(core.ids.AtomId.fromIndex(1), constraints[0].double_b);
}

test "a small ring contributes nothing even with a stereo bond" {
    const ring_size = 6;
    var atoms: [ring_size]model.Atom = undefined;
    var bonds: [ring_size]model.Bond = undefined;
    for (&atoms, &bonds, 0..) |*atom, *bond, index| {
        atom.* = .{
            .id = core.ids.AtomId.fromIndex(@intCast(index)),
            .input_index = @intCast(index),
            .atomic_number = .carbon,
        };
        bond.* = .{
            .id = core.ids.BondId.fromIndex(@intCast(index)),
            .input_index = @intCast(index),
            .start = core.ids.AtomId.fromIndex(@intCast(index)),
            .end = core.ids.AtomId.fromIndex(@intCast((index + 1) % ring_size)),
            .input_order = .single,
            .effective_order = .single,
        };
    }
    bonds[0].effective_order = .double;
    bonds[0].stereo = .cis;

    var graph = try topology.Graph.init(testing.allocator, &atoms, &bonds);
    defer graph.deinit();
    var membership = try topology.RingMembership.init(testing.allocator, graph, &bonds);
    defer membership.deinit();

    const constraints = try buildChiralInversionConstraints(
        testing.allocator,
        &atoms,
        &bonds,
        graph,
        membership,
    );
    defer testing.allocator.free(constraints);
    try testing.expectEqual(@as(usize, 0), constraints.len);
}

/// Peptide bond inversion constraints for the DISCRETE pass
/// (`addPeptideBondInversionConstraintsOfMolecule`). Upstream adds these to
/// the interaction set that `scoreClashes` and therefore `flipFragments` read,
/// not to the minimization set, so they influence the clean-pose verdict.
///
/// Upstream walks `std::set<sketcherMinimizerAtom*>`, whose iteration order is
/// pointer order and therefore allocation-address order. This walks atoms in
/// index order instead. That is a correction rather than an emulation, in the
/// same family as the other heap-order corrections this port makes: the
/// emitted set is identical, only its order is made deterministic.
pub fn buildPeptideBondInversionConstraints(
    allocator: std.mem.Allocator,
    atoms: []const model.Atom,
    bonds: []const model.Bond,
    graph: topology.Graph,
) core.errors.Error![]core.interaction.EzConstraint {
    const cheto_carbons = allocator.alloc(bool, atoms.len) catch return error.OutOfMemory;
    defer allocator.free(cheto_carbons);
    const amino_nitrogens = allocator.alloc(bool, atoms.len) catch return error.OutOfMemory;
    defer allocator.free(amino_nitrogens);
    const alpha_carbons = allocator.alloc(bool, atoms.len) catch return error.OutOfMemory;
    defer allocator.free(alpha_carbons);
    @memset(cheto_carbons, false);
    @memset(amino_nitrogens, false);
    @memset(alpha_carbons, false);

    var cheto_count: usize = 0;
    var amino_count: usize = 0;
    var alpha_count: usize = 0;

    for (atoms) |atom| {
        const index = atom.id.index();
        if (atom.atomic_number == .nitrogen) {
            amino_nitrogens[index] = true;
            amino_count += 1;
            continue;
        }
        if (atom.atomic_number != .carbon) continue;
        if (try isChetoCarbon(atom.id, atoms, bonds, graph)) {
            cheto_carbons[index] = true;
            cheto_count += 1;
        }
    }
    if (cheto_count < 2 or amino_count < 2) return allocator.alloc(core.interaction.EzConstraint, 0) catch
        return error.OutOfMemory;

    for (atoms) |atom| {
        const index = atom.id.index();
        if (atom.atomic_number != .carbon or cheto_carbons[index]) continue;
        var bonded_to_cheto = false;
        var bonded_to_amino = false;
        for (graph.neighbors(atom.id)) |neighbor| {
            if (cheto_carbons[neighbor.index()]) bonded_to_cheto = true;
            if (amino_nitrogens[neighbor.index()]) bonded_to_amino = true;
        }
        if (bonded_to_cheto and bonded_to_amino) {
            alpha_carbons[index] = true;
            alpha_count += 1;
        }
    }
    if (alpha_count < 2) return allocator.alloc(core.interaction.EzConstraint, 0) catch
        return error.OutOfMemory;

    var constraints: std.ArrayList(core.interaction.EzConstraint) = .empty;
    defer constraints.deinit(allocator);
    const sequences = [_][4][]const bool{
        .{ cheto_carbons, amino_nitrogens, alpha_carbons, cheto_carbons },
        .{ amino_nitrogens, alpha_carbons, cheto_carbons, amino_nitrogens },
        .{ alpha_carbons, cheto_carbons, amino_nitrogens, alpha_carbons },
    };
    for (sequences) |sequence| {
        try appendSequenceMatches(allocator, &constraints, atoms, graph, sequence);
    }
    return constraints.toOwnedSlice(allocator) catch return error.OutOfMemory;
}

/// Four bonded atoms drawn from the four sets in order, with no distinctness
/// requirement - upstream imposes none.
fn appendSequenceMatches(
    allocator: std.mem.Allocator,
    constraints: *std.ArrayList(core.interaction.EzConstraint),
    atoms: []const model.Atom,
    graph: topology.Graph,
    sets: [4][]const bool,
) core.errors.Error!void {
    for (atoms) |atom| {
        if (!sets[0][atom.id.index()]) continue;
        for (graph.neighbors(atom.id)) |second| {
            if (!sets[1][second.index()]) continue;
            for (graph.neighbors(second)) |third| {
                if (!sets[2][third.index()]) continue;
                for (graph.neighbors(third)) |fourth| {
                    if (!sets[3][fourth.index()]) continue;
                    constraints.append(allocator, .{
                        .side_a = atom.id,
                        .double_a = second,
                        .double_b = third,
                        .side_b = fourth,
                        // Upstream passes cis = false for every peptide group.
                        .is_z = false,
                    }) catch return error.OutOfMemory;
                }
            }
        }
    }
}

test "a dipeptide backbone yields peptide constraints and a plain chain yields none" {
    // C(=O)-N-C(alpha)-C(=O)-N-C(alpha): two cheto carbons, two nitrogens,
    // two alpha carbons.
    const atoms = [_]model.Atom{
        .{ .id = core.ids.AtomId.fromIndex(0), .input_index = 0, .atomic_number = .carbon },
        .{ .id = core.ids.AtomId.fromIndex(1), .input_index = 1, .atomic_number = .oxygen },
        .{ .id = core.ids.AtomId.fromIndex(2), .input_index = 2, .atomic_number = .nitrogen },
        .{ .id = core.ids.AtomId.fromIndex(3), .input_index = 3, .atomic_number = .carbon },
        .{ .id = core.ids.AtomId.fromIndex(4), .input_index = 4, .atomic_number = .carbon },
        .{ .id = core.ids.AtomId.fromIndex(5), .input_index = 5, .atomic_number = .oxygen },
        .{ .id = core.ids.AtomId.fromIndex(6), .input_index = 6, .atomic_number = .nitrogen },
        .{ .id = core.ids.AtomId.fromIndex(7), .input_index = 7, .atomic_number = .carbon },
        .{ .id = core.ids.AtomId.fromIndex(8), .input_index = 8, .atomic_number = .carbon },
        .{ .id = core.ids.AtomId.fromIndex(9), .input_index = 9, .atomic_number = .oxygen },
    };
    const bonds = [_]model.Bond{
        peptideBond(0, 0, 1, .double),
        peptideBond(1, 0, 2, .single),
        peptideBond(2, 2, 3, .single),
        peptideBond(3, 3, 4, .single),
        peptideBond(4, 4, 5, .double),
        peptideBond(5, 4, 6, .single),
        peptideBond(6, 6, 7, .single),
        peptideBond(7, 7, 8, .single),
        peptideBond(8, 8, 9, .double),
    };
    var graph = try topology.Graph.init(testing.allocator, &atoms, &bonds);
    defer graph.deinit();

    const constraints = try buildPeptideBondInversionConstraints(
        testing.allocator,
        &atoms,
        &bonds,
        graph,
    );
    defer testing.allocator.free(constraints);
    try testing.expect(constraints.len > 0);
    for (constraints) |constraint| try testing.expect(!constraint.is_z);

    // The same skeleton without nitrogens has no peptide backbone at all.
    var carbons = atoms;
    carbons[2].atomic_number = .carbon;
    carbons[6].atomic_number = .carbon;
    var plain_graph = try topology.Graph.init(testing.allocator, &carbons, &bonds);
    defer plain_graph.deinit();
    const none = try buildPeptideBondInversionConstraints(
        testing.allocator,
        &carbons,
        &bonds,
        plain_graph,
    );
    defer testing.allocator.free(none);
    try testing.expectEqual(@as(usize, 0), none.len);
}

fn peptideBond(index: u32, start: u32, end: u32, order: core.chemistry.BondOrder) model.Bond {
    return .{
        .id = core.ids.BondId.fromIndex(index),
        .input_index = index,
        .start = core.ids.AtomId.fromIndex(start),
        .end = core.ids.AtomId.fromIndex(end),
        .input_order = order,
        .effective_order = order,
    };
}
