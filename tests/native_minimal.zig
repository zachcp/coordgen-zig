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
    for (first.coordinates) |coordinate| try std.testing.expect(coordinate.isFinite());
    for (bonds) |bond| {
        const delta_x = first.coordinates[bond.start].x - first.coordinates[bond.end].x;
        const delta_y = first.coordinates[bond.start].y - first.coordinates[bond.end].y;
        try std.testing.expectApproxEqAbs(api.bond_length, @sqrt(delta_x * delta_x + delta_y * delta_y), 0.001);
    }
    // All-zero preparation output (layout omitted) fails this executable fact.
    try std.testing.expect(first.coordinates[0].x != first.coordinates[1].x or first.coordinates[0].y != first.coordinates[1].y);
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
    // Basic layout anchors canonical atom zero at the origin. This assertion
    // fails if internal coordinates are copied directly or either map is used
    // in the transposed direction.
    try std.testing.expectEqual(api.Vec2{}, result.coordinates[result.internal_to_input[0]]);
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

test "minimal native generation cleans every injected allocation failure" {
    const atoms = [_]api.AtomInput{ .{}, .{}, .{} };
    const bonds = [_]api.BondInput{ .{ .start = 0, .end = 1 }, .{ .start = 1, .end = 2 } };
    try std.testing.checkAllAllocationFailures(
        std.testing.allocator,
        generateAndDiscard,
        .{api.Input{ .atoms = &atoms, .bonds = &bonds }},
    );
}

test "minimal native generation explicitly rejects domains owned by later phases" {
    const atoms = [_]api.AtomInput{ .{}, .{}, .{} };
    const path = [_]api.BondInput{ .{ .start = 0, .end = 1 }, .{ .start = 1, .end = 2 } };
    const cycle = [_]api.BondInput{ path[0], path[1], .{ .start = 2, .end = 0 } };
    const disconnected = [_]api.BondInput{path[0]};
    var cycle_result = try generate(std.testing.allocator, .{ .atoms = &atoms, .bonds = &cycle });
    cycle_result.deinit();
    try expectUnsupported(.{ .atoms = &atoms, .bonds = &disconnected });

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
    changed_bonds[0].order = .zero;
    try expectUnsupported(.{ .atoms = &atoms, .bonds = &changed_bonds });
    changed_bonds = path;
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
        .residues = &.{.{ .atom = 0 }},
    });
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
        api.Options{ .build_from_fragments = true },
        api.Options{ .debug_coordinates = true },
        api.Options{ .template_directory = "fixtures" },
    }) |options| try expectUnsupported(.{ .atoms = &atoms, .bonds = &path, .options = options });
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
    try std.testing.checkAllAllocationFailures(std.testing.allocator, generateAndDiscard, .{input});
    try std.testing.checkAllAllocationFailures(std.testing.allocator, generateAndDiscard, .{api.Input{
        .atoms = &atoms,
        .bonds = &bonds,
        .options = .{ .force_open_macrocycles = true },
    }});
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
    try std.testing.checkAllAllocationFailures(std.testing.allocator, generateAndDiscard, .{input});
}

test "minimal native validation rejects malformed input before generation" {
    const atoms = [_]api.AtomInput{ .{}, .{} };
    try std.testing.expectError(error.InvalidAtomIndex, generate(std.testing.allocator, .{
        .atoms = &atoms,
        .bonds = &.{.{ .start = 0, .end = 2 }},
    }));
}
