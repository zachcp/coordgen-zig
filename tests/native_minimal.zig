const std = @import("std");
const api = @import("api");
const minimal = @import("generator");

fn generate(allocator: std.mem.Allocator, input: api.Input) !minimal.Result {
    try input.validate();
    return minimal.generateValidated(allocator, input);
}

test "minimal native ethane and propane are finite deterministic caller-order layouts" {
    const atoms = [_]api.AtomInput{ .{}, .{}, .{} };
    const bonds = [_]api.BondInput{ .{ .start = 0, .end = 1 }, .{ .start = 1, .end = 2 } };
    const input = api.Input{ .atoms = &atoms, .bonds = &bonds };
    var first = try generate(std.testing.allocator, input);
    defer first.deinit();
    var second = try generate(std.testing.allocator, input);
    defer second.deinit();
    try std.testing.expectEqualSlices(api.Vec2, first.coordinates, second.coordinates);
    try std.testing.expect(first.clean_pose);
    for (first.coordinates) |coordinate| try std.testing.expect(coordinate.isFinite());
    for (bonds) |bond| {
        const delta_x = first.coordinates[bond.start].x - first.coordinates[bond.end].x;
        const delta_y = first.coordinates[bond.start].y - first.coordinates[bond.end].y;
        try std.testing.expectApproxEqAbs(api.bond_length, @sqrt(delta_x * delta_x + delta_y * delta_y), 0.001);
    }
    // All-zero preparation output (layout omitted) fails this executable fact.
    try std.testing.expect(first.coordinates[0].x != first.coordinates[1].x or first.coordinates[0].y != first.coordinates[1].y);
}

test "five-atom path clean pose matches the stable upstream probe" {
    const atoms = [_]api.AtomInput{ .{}, .{}, .{}, .{}, .{} };
    const bonds = [_]api.BondInput{
        .{ .start = 0, .end = 1 },
        .{ .start = 1, .end = 2 },
        .{ .start = 2, .end = 3 },
        .{ .start = 3, .end = 4 },
    };
    var result = try generate(std.testing.allocator, .{ .atoms = &atoms, .bonds = &bonds });
    defer result.deinit();
    try std.testing.expect(!result.clean_pose);
}

test "minimal native result maps canonical coordinates back to caller order without mutating input" {
    const atoms = [_]api.AtomInput{
        .{ .atomic_number = .carbon },
        .{ .atomic_number = .oxygen },
        .{ .atomic_number = .nitrogen },
    };
    const bonds = [_]api.BondInput{ .{ .start = 0, .end = 1 }, .{ .start = 1, .end = 2 } };
    const atoms_before = atoms;
    const bonds_before = bonds;
    var result = try generate(std.testing.allocator, .{ .atoms = &atoms, .bonds = &bonds });
    defer result.deinit();
    try std.testing.expectEqualDeep(atoms_before, atoms);
    try std.testing.expectEqualDeep(bonds_before, bonds);
    try std.testing.expect(!std.mem.eql(u32, result.input_to_internal, &.{ 0, 1, 2 }));
    for (result.internal_to_input, 0..) |input_index, internal_index| {
        try std.testing.expectEqual(@as(u32, @intCast(internal_index)), result.input_to_internal[input_index]);
    }
    for (bonds) |bond| {
        const delta_x = result.coordinates[bond.start].x - result.coordinates[bond.end].x;
        const delta_y = result.coordinates[bond.start].y - result.coordinates[bond.end].y;
        try std.testing.expectApproxEqAbs(api.bond_length, @sqrt(delta_x * delta_x + delta_y * delta_y), 0.001);
    }
}

fn generateAndDiscard(allocator: std.mem.Allocator, input: api.Input) !void {
    var result = try generate(allocator, input);
    result.deinit();
}

fn expectUnsupported(input: api.Input) !void {
    try std.testing.expectError(error.Unsupported, generate(std.testing.allocator, input));
}

/// Fail every allocation of one generation in turn, and prove the count is
/// stable rather than assuming it.
///
/// Every measurement runs against a FRESH arena. Counting against a shared
/// allocator is not reproducible: whether an ArrayList growth is served by an
/// in-place resize or by a new allocation depends on residual heap layout from
/// whatever ran before, so the same input can measure 108 allocations on one
/// run and 109 on the next with no change in behaviour. That oscillation is a
/// property of the measurement, not of the code, and it is what made this
/// helper necessary in the first place.
///
/// `backing` is still swept for leaks: each run asserts allocated == freed,
/// and the caller's std.testing.allocator would fail the binary on any leak
/// the arena hid.
fn checkStableProteinAllocations(backing: std.mem.Allocator, input: api.Input) !void {
    var leak_check = std.testing.FailingAllocator.init(backing, .{});
    try generateAndDiscard(leak_check.allocator(), input);
    try std.testing.expectEqual(leak_check.allocated_bytes, leak_check.freed_bytes);

    var counts: [2]usize = undefined;
    for (&counts) |*count| {
        var arena = std.heap.ArenaAllocator.init(std.heap.page_allocator);
        defer arena.deinit();
        var measured = std.testing.FailingAllocator.init(arena.allocator(), .{});
        try generateAndDiscard(measured.allocator(), input);
        count.* = measured.alloc_index;
    }
    try std.testing.expectEqual(counts[0], counts[1]);

    for (0..counts[1]) |index| {
        var arena = std.heap.ArenaAllocator.init(std.heap.page_allocator);
        defer arena.deinit();
        var failing = std.testing.FailingAllocator.init(arena.allocator(), .{ .fail_index = index });
        if (generateAndDiscard(failing.allocator(), input)) |_| {
            if (failing.has_induced_failure) return error.SwallowedOutOfMemoryError;
            return error.NondeterministicMemoryUsage;
        } else |err| {
            if (err != error.OutOfMemory) return err;
            if (failing.allocated_bytes != failing.freed_bytes) return error.MemoryLeakDetected;
        }
    }
}

test "minimal native generation cleans every injected allocation failure" {
    const atoms = [_]api.AtomInput{ .{}, .{}, .{} };
    const bonds = [_]api.BondInput{ .{ .start = 0, .end = 1 }, .{ .start = 1, .end = 2 } };
    // Same coverage as checkAllAllocationFailures - every allocation index is
    // failed in turn - with the ReleaseSafe warm-up allocation discharged
    // first and stability proved rather than assumed. cgz-7v2.21 moved the
    // result allocations ahead of the pipeline so both entry points share one
    // implementation, which brought this path into the same wrapper warm-up
    // the residue paths already document.
    try checkStableProteinAllocations(
        std.testing.allocator,
        .{ .atoms = &atoms, .bonds = &bonds },
    );
}

test "minimal native generation explicitly rejects domains owned by later phases" {
    const atoms = [_]api.AtomInput{ .{}, .{}, .{} };
    const path = [_]api.BondInput{ .{ .start = 0, .end = 1 }, .{ .start = 1, .end = 2 } };
    const cycle = [_]api.BondInput{ path[0], path[1], .{ .start = 2, .end = 0 } };
    var cycle_result = try generate(std.testing.allocator, .{ .atoms = &atoms, .bonds = &cycle });
    cycle_result.deinit();

    var changed = atoms;
    changed[0].hidden = true;
    try expectUnsupported(.{ .atoms = &changed, .bonds = &path });
    changed = atoms;
    changed[0].fixed = true;
    try expectUnsupported(.{ .atoms = &changed, .bonds = &path });
    changed = atoms;
    changed[0].constrained = true;
    try expectUnsupported(.{ .atoms = &changed, .bonds = &path });
    changed = atoms;
    changed[0].template_coordinates = .{};
    try expectUnsupported(.{ .atoms = &changed, .bonds = &path });
    changed = atoms;
    changed[0].coordinates_3d = .{};
    try expectUnsupported(.{ .atoms = &changed, .bonds = &path });
    changed = atoms;
    changed[0].stereo.value = .r;
    try expectUnsupported(.{ .atoms = &changed, .bonds = &path });

    var changed_bonds = path;
    changed_bonds[0].skip = true;
    try expectUnsupported(.{ .atoms = &atoms, .bonds = &changed_bonds });
    changed_bonds = path;
    changed_bonds[0].stereo.value = .z;
    try expectUnsupported(.{ .atoms = &atoms, .bonds = &changed_bonds });
    changed_bonds = path;
    changed_bonds[0].display = .solid_forward;
    try expectUnsupported(.{ .atoms = &atoms, .bonds = &changed_bonds });

    try expectUnsupported(.{
        .atoms = &atoms,
        .bonds = &path,
        .residue_interactions = &.{.{ .start = 0, .end = 2 }},
    });
    try expectUnsupported(.{
        .atoms = &atoms,
        .bonds = &path,
        .extra_bonds = &.{.{ .start = 0, .end = 2 }},
    });

    inline for (.{
        api.Options{ .even_angles = true },
        api.Options{ .constrain_all_atoms = true },
        api.Options{ .debug_coordinates = true },
        api.Options{ .template_directory = "fixtures" },
    }) |options| try expectUnsupported(.{ .atoms = &atoms, .bonds = &path, .options = options });
}

test "minimal native generation places residue representatives around a ligand" {
    const atoms = [_]api.AtomInput{ .{}, .{}, .{}, .{} };
    const bonds = [_]api.BondInput{ .{ .start = 0, .end = 1 }, .{ .start = 1, .end = 2 } };
    const residues = [_]api.ResidueInput{.{
        .atom = 3,
        .chain = "A",
        .residue_number = 7,
        .closest_ligand_atom = 1,
    }};
    const interactions = [_]api.ResidueInteractionInput{.{ .start = 3, .end = 1 }};
    const input = api.Input{
        .atoms = &atoms,
        .bonds = &bonds,
        .residues = &residues,
        .residue_interactions = &interactions,
    };
    var first = try generate(std.testing.allocator, input);
    defer first.deinit();
    var repeated = try generate(std.testing.allocator, input);
    defer repeated.deinit();
    try std.testing.expectEqualSlices(api.Vec2, first.coordinates, repeated.coordinates);
    const dx = first.coordinates[3].x - first.coordinates[1].x;
    const dy = first.coordinates[3].y - first.coordinates[1].y;
    try std.testing.expect(@sqrt(dx * dx + dy * dy) > 40);
    for (first.coordinates) |coordinate| try std.testing.expect(coordinate.isFinite());
    try checkStableProteinAllocations(std.testing.allocator, input);
}

test "minimal native generation places interacting protein-only chains" {
    const atoms = [_]api.AtomInput{ .{}, .{}, .{}, .{} };
    const residues = [_]api.ResidueInput{
        .{ .atom = 0, .chain = "A", .residue_number = 1 },
        .{ .atom = 1, .chain = "A", .residue_number = 2 },
        .{ .atom = 2, .chain = "B", .residue_number = 1 },
        .{ .atom = 3, .chain = "B", .residue_number = 2 },
    };
    const interactions = [_]api.ResidueInteractionInput{
        .{ .start = 0, .end = 2 },
        .{ .start = 1, .end = 3 },
    };
    const input = api.Input{
        .atoms = &atoms,
        .bonds = &.{},
        .residues = &residues,
        .residue_interactions = &interactions,
    };
    var first = try generate(std.testing.allocator, input);
    defer first.deinit();
    var repeated = try generate(std.testing.allocator, input);
    defer repeated.deinit();
    try std.testing.expectEqualSlices(api.Vec2, first.coordinates, repeated.coordinates);
    for (first.coordinates) |coordinate| try std.testing.expect(coordinate.isFinite());
    try std.testing.expect(!std.meta.eql(first.coordinates[0], first.coordinates[1]));
    try std.testing.expect(!std.meta.eql(first.coordinates[0], first.coordinates[2]));
    try checkStableProteinAllocations(std.testing.allocator, input);
}

test "minimal native generation composes ligand proximity and residue placement" {
    const atoms = [_]api.AtomInput{ .{}, .{}, .{}, .{} };
    const bonds = [_]api.BondInput{
        .{ .start = 0, .end = 1 },
        .{ .start = 1, .end = 2, .order = .zero },
    };
    const residues = [_]api.ResidueInput{.{
        .atom = 3,
        .chain = "A",
        .residue_number = 4,
        .closest_ligand_atom = 1,
    }};
    const interactions = [_]api.ResidueInteractionInput{.{ .start = 3, .end = 1 }};
    const input = api.Input{
        .atoms = &atoms,
        .bonds = &bonds,
        .residues = &residues,
        .residue_interactions = &interactions,
    };
    var first = try generate(std.testing.allocator, input);
    defer first.deinit();
    var repeated = try generate(std.testing.allocator, input);
    defer repeated.deinit();
    try std.testing.expectEqualSlices(api.Vec2, first.coordinates, repeated.coordinates);
    for (first.coordinates) |coordinate| try std.testing.expect(coordinate.isFinite());
    try std.testing.expect(!std.meta.eql(first.coordinates[1], first.coordinates[2]));
    try std.testing.expect(!std.meta.eql(first.coordinates[1], first.coordinates[3]));
}

test "residue generation is invariant under reversed caller order" {
    const atoms = [_]api.AtomInput{
        .{ .atomic_number = .carbon },
        .{ .atomic_number = .oxygen },
        .{ .atomic_number = .nitrogen },
        .{ .atomic_number = .sulfur },
    };
    const bonds = [_]api.BondInput{ .{ .start = 0, .end = 1 }, .{ .start = 1, .end = 2 } };
    const residues = [_]api.ResidueInput{.{ .atom = 3, .chain = "A", .residue_number = 9, .closest_ligand_atom = 1 }};
    const interactions = [_]api.ResidueInteractionInput{.{ .start = 3, .end = 1 }};
    var forward = try generate(std.testing.allocator, .{
        .atoms = &atoms,
        .bonds = &bonds,
        .residues = &residues,
        .residue_interactions = &interactions,
    });
    defer forward.deinit();

    const reversed_atoms = [_]api.AtomInput{ atoms[3], atoms[2], atoms[1], atoms[0] };
    const reversed_bonds = [_]api.BondInput{ .{ .start = 3, .end = 2 }, .{ .start = 2, .end = 1 } };
    const reversed_residues = [_]api.ResidueInput{.{ .atom = 0, .chain = "A", .residue_number = 9, .closest_ligand_atom = 2 }};
    const reversed_interactions = [_]api.ResidueInteractionInput{.{ .start = 0, .end = 2 }};
    var reversed = try generate(std.testing.allocator, .{
        .atoms = &reversed_atoms,
        .bonds = &reversed_bonds,
        .residues = &reversed_residues,
        .residue_interactions = &reversed_interactions,
    });
    defer reversed.deinit();
    for (forward.coordinates, 0..) |coordinate, index| {
        try std.testing.expectEqual(coordinate, reversed.coordinates[forward.coordinates.len - 1 - index]);
    }
    try std.testing.expectEqual(forward.clean_pose, reversed.clean_pose);
}

const ThreadGenerationContext = struct {
    input: api.Input,
    output: [4]api.Vec2 = undefined,
    clean_pose: bool = false,
    failure: ?anyerror = null,
};

fn generateInThread(context: *ThreadGenerationContext) void {
    var result = generate(std.heap.page_allocator, context.input) catch |err| {
        context.failure = err;
        return;
    };
    defer result.deinit();
    @memcpy(&context.output, result.coordinates);
    context.clean_pose = result.clean_pose;
}

test "residue generation has independent simultaneous contexts" {
    const atoms = [_]api.AtomInput{ .{}, .{}, .{}, .{} };
    const bonds = [_]api.BondInput{ .{ .start = 0, .end = 1 }, .{ .start = 1, .end = 2 } };
    const residues = [_]api.ResidueInput{.{ .atom = 3, .chain = "A", .residue_number = 2, .closest_ligand_atom = 1 }};
    const interactions = [_]api.ResidueInteractionInput{.{ .start = 3, .end = 1 }};
    const input = api.Input{
        .atoms = &atoms,
        .bonds = &bonds,
        .residues = &residues,
        .residue_interactions = &interactions,
    };
    var first = ThreadGenerationContext{ .input = input };
    var second = ThreadGenerationContext{ .input = input };
    const first_thread = try std.Thread.spawn(.{}, generateInThread, .{&first});
    const second_thread = try std.Thread.spawn(.{}, generateInThread, .{&second});
    first_thread.join();
    second_thread.join();
    if (first.failure) |err| return err;
    if (second.failure) |err| return err;
    try std.testing.expectEqual(first.output, second.output);
    try std.testing.expectEqual(first.clean_pose, second.clean_pose);
}

test "minimal native generation reaches macrocycle lattice and forced-open fallback" {
    const ring_size = 13;
    var atoms: [ring_size]api.AtomInput = undefined;
    var bonds: [ring_size]api.BondInput = undefined;
    for (&atoms, &bonds, 0..) |*atom, *bond, index| {
        atom.* = .{};
        bond.* = .{
            .start = @intCast(index),
            .end = @intCast((index + 1) % ring_size),
        };
    }
    const input = api.Input{ .atoms = &atoms, .bonds = &bonds };
    var lattice = try generate(std.testing.allocator, input);
    defer lattice.deinit();
    var repeated = try generate(std.testing.allocator, input);
    defer repeated.deinit();
    try std.testing.expectEqualSlices(api.Vec2, lattice.coordinates, repeated.coordinates);
    for (lattice.coordinates) |coordinate| try std.testing.expect(coordinate.isFinite());

    var opened = try generate(std.testing.allocator, .{
        .atoms = &atoms,
        .bonds = &bonds,
        .options = .{ .force_open_macrocycles = true },
    });
    defer opened.deinit();
    var changed = false;
    for (lattice.coordinates, opened.coordinates) |lattice_coordinate, opened_coordinate| {
        changed = changed or !std.meta.eql(lattice_coordinate, opened_coordinate);
    }
    try std.testing.expect(changed);
    for (opened.coordinates) |coordinate| try std.testing.expect(coordinate.isFinite());
    try checkStableProteinAllocations(std.testing.allocator, input);
    try checkStableProteinAllocations(std.testing.allocator, .{
        .atoms = &atoms,
        .bonds = &bonds,
        .options = .{ .force_open_macrocycles = true },
    });
}

test "minimal native generation runs discrete search for macrocycle substituents" {
    const ring_size = 13;
    var atoms: [ring_size + 2]api.AtomInput = undefined;
    @memset(&atoms, .{});
    var bonds: [ring_size + 2]api.BondInput = undefined;
    for (bonds[0..ring_size], 0..) |*bond, index| bond.* = .{
        .start = @intCast(index),
        .end = @intCast((index + 1) % ring_size),
    };
    bonds[ring_size] = .{ .start = 2, .end = ring_size };
    bonds[ring_size + 1] = .{ .start = ring_size, .end = ring_size + 1 };
    const input = api.Input{ .atoms = &atoms, .bonds = &bonds };
    var first = try generate(std.testing.allocator, input);
    defer first.deinit();
    var second = try generate(std.testing.allocator, input);
    defer second.deinit();
    try std.testing.expectEqualSlices(api.Vec2, first.coordinates, second.coordinates);
    for (first.coordinates) |coordinate| try std.testing.expect(coordinate.isFinite());
    for (bonds[ring_size..]) |bond| {
        const dx = first.coordinates[bond.start].x - first.coordinates[bond.end].x;
        const dy = first.coordinates[bond.start].y - first.coordinates[bond.end].y;
        try std.testing.expectApproxEqAbs(api.bond_length, @sqrt(dx * dx + dy * dy), 0.01);
    }
    try checkStableProteinAllocations(std.testing.allocator, input);
}

test "minimal native generation arranges disconnected neutral components" {
    const atoms = [_]api.AtomInput{ .{}, .{}, .{}, .{} };
    const bonds = [_]api.BondInput{
        .{ .start = 0, .end = 1 },
        .{ .start = 2, .end = 3 },
    };
    const input = api.Input{ .atoms = &atoms, .bonds = &bonds };
    var result = try generate(std.testing.allocator, input);
    defer result.deinit();
    for (bonds) |bond| {
        const dx = result.coordinates[bond.start].x - result.coordinates[bond.end].x;
        const dy = result.coordinates[bond.start].y - result.coordinates[bond.end].y;
        try std.testing.expectApproxEqAbs(api.bond_length, @sqrt(dx * dx + dy * dy), 0.001);
    }
    for (result.coordinates[0..2]) |first| {
        for (result.coordinates[2..4]) |second| {
            try std.testing.expect(@abs(first.x - second.x) >= api.bond_length or
                @abs(first.y - second.y) >= api.bond_length);
        }
    }
    try std.testing.checkAllAllocationFailures(std.testing.allocator, generateAndDiscard, .{input});
}

test "minimal native generation places an acyclic proximity child from a large center" {
    var atoms: [10]api.AtomInput = undefined;
    @memset(&atoms, .{});
    var bonds: [9]api.BondInput = undefined;
    for (bonds[0..7], 0..) |*bond, index| bond.* = .{
        .start = @intCast(index),
        .end = @intCast(index + 1),
    };
    bonds[7] = .{ .start = 8, .end = 9 };
    bonds[8] = .{ .start = 7, .end = 8, .order = .zero };
    const input = api.Input{ .atoms = &atoms, .bonds = &bonds };
    var result = try generate(std.testing.allocator, input);
    defer result.deinit();
    const dx = result.coordinates[7].x - result.coordinates[8].x;
    const dy = result.coordinates[7].y - result.coordinates[8].y;
    try std.testing.expect(@sqrt(dx * dx + dy * dy) >= api.bond_length * 2);
    try std.testing.checkAllAllocationFailures(std.testing.allocator, generateAndDiscard, .{input});
}

test "minimal native generation uses general placement for a small proximity pair" {
    const atoms = [_]api.AtomInput{ .{}, .{}, .{} };
    const bonds = [_]api.BondInput{
        .{ .start = 0, .end = 1 },
        .{ .start = 1, .end = 2, .order = .zero },
    };
    const input = api.Input{ .atoms = &atoms, .bonds = &bonds };
    var result = try generate(std.testing.allocator, input);
    defer result.deinit();
    for (result.coordinates) |coordinate| try std.testing.expect(coordinate.isFinite());
    const dx = result.coordinates[1].x - result.coordinates[2].x;
    const dy = result.coordinates[1].y - result.coordinates[2].y;
    try std.testing.expect(@sqrt(dx * dx + dy * dy) >= api.bond_length - 0.001);
}

test "minimal native generation lays out acyclic and cyclic proximity meta graphs" {
    const atoms = [_]api.AtomInput{ .{}, .{}, .{} };
    const path_bonds = [_]api.BondInput{
        .{ .start = 0, .end = 1, .order = .zero },
        .{ .start = 1, .end = 2, .order = .zero },
    };
    var path = try generate(std.testing.allocator, .{ .atoms = &atoms, .bonds = &path_bonds });
    defer path.deinit();
    for (path.coordinates) |coordinate| try std.testing.expect(coordinate.isFinite());
    try std.testing.expect(!std.meta.eql(path.coordinates[0], path.coordinates[2]));

    const cycle_bonds = [_]api.BondInput{
        path_bonds[0],
        path_bonds[1],
        .{ .start = 2, .end = 0, .order = .zero },
    };
    var cycle = try generate(std.testing.allocator, .{ .atoms = &atoms, .bonds = &cycle_bonds });
    defer cycle.deinit();
    for (cycle.coordinates) |coordinate| try std.testing.expect(coordinate.isFinite());
    for (cycle_bonds) |bond| {
        const dx = cycle.coordinates[bond.start].x - cycle.coordinates[bond.end].x;
        const dy = cycle.coordinates[bond.start].y - cycle.coordinates[bond.end].y;
        try std.testing.expectApproxEqAbs(api.bond_length, @sqrt(dx * dx + dy * dy), 0.001);
    }
}

test "minimal native validation rejects malformed input before generation" {
    const atoms = [_]api.AtomInput{ .{}, .{} };
    try std.testing.expectError(error.InvalidAtomIndex, generate(std.testing.allocator, .{
        .atoms = &atoms,
        .bonds = &.{.{ .start = 0, .end = 2 }},
    }));
}

