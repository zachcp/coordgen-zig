const std = @import("std");
const api = @import("api");
const conformance = @import("conformance");
const generator = @import("generator");

pub fn main(init: std.process.Init) !void {
    const args = try init.minimal.args.toSlice(init.arena.allocator());
    if (args.len != 2) return error.InvalidArguments;
    const source = try std.Io.Dir.cwd().readFileAlloc(
        init.io,
        args[1],
        init.gpa,
        .limited(1 << 20),
    );
    defer init.gpa.free(source);
    var reading = try conformance.mae.parse(init.gpa, source, null);
    defer reading.deinit();
    if (reading.structures.len != 1) return error.InvalidFixture;
    const structure = reading.structures[0];
    const atoms = try init.gpa.alloc(api.AtomInput, structure.atoms.len);
    defer init.gpa.free(atoms);
    for (atoms, structure.atoms) |*atom, source_atom| atom.* = .{
        .atomic_number = api.AtomicNumber.fromPublic(source_atom.atomic_number) orelse return error.InvalidFixture,
    };
    const bonds = try init.gpa.alloc(api.BondInput, structure.bonds.len);
    defer init.gpa.free(bonds);
    for (bonds, structure.bonds) |*bond, source_bond| bond.* = .{
        .start = source_bond.from,
        .end = source_bond.to,
        .order = api.BondOrder.fromInt(source_bond.order) orelse return error.InvalidFixture,
    };
    const input = api.Input{ .atoms = atoms, .bonds = bonds };
    try input.validate();
    var first = try generator.generateValidated(init.gpa, input);
    defer first.deinit();
    var repeated = try generator.generateValidated(init.gpa, input);
    defer repeated.deinit();
    if (first.coordinates.len != structure.atoms.len or repeated.coordinates.len != first.coordinates.len) return error.InvalidResult;
    for (first.coordinates, repeated.coordinates) |coordinate, repeated_coordinate| {
        if (!coordinate.isFinite() or !std.meta.eql(coordinate, repeated_coordinate)) return error.InvalidResult;
    }
    var opened = try generator.generateValidated(init.gpa, api.Input{
        .atoms = atoms,
        .bonds = bonds,
        .options = .{ .force_open_macrocycles = true },
    });
    defer opened.deinit();
    for (opened.coordinates) |coordinate| if (!coordinate.isFinite()) return error.InvalidResult;
}
