const std = @import("std");
const coordgen = @import("coordgen");

pub fn main() !void {
    var gpa: std.heap.DebugAllocator(.{}) = .init;
    defer if (gpa.deinit() == .leak) @panic("installed Zig consumer leaked its result");

    const atoms = [_]coordgen.AtomInput{
        .{},
        .{ .atomic_number = .oxygen },
    };
    const bonds = [_]coordgen.BondInput{
        .{ .start = 0, .end = 1, .order = .double },
    };
    const input: coordgen.Input = .{ .atoms = &atoms, .bonds = &bonds };
    try input.validate();
    var result = try coordgen.api.generate(gpa.allocator(), input);
    defer result.deinit();
    if (result.coordinates.len != atoms.len or
        result.effective_bond_orders.len != bonds.len or
        !result.coordinates[0].isFinite() or
        !result.coordinates[1].isFinite())
    {
        return error.InvalidGeneratedResult;
    }

    if (input.options.precision != coordgen.api.Precision.standard) {
        return error.UnexpectedDefaultPrecision;
    }
    if (coordgen.bond_length != 50.0) return error.UnexpectedBondLength;

    const invalid_bonds = [_]coordgen.BondInput{
        .{ .start = 0, .end = 2 },
    };
    const invalid: coordgen.Input = .{ .atoms = &atoms, .bonds = &invalid_bonds };
    invalid.validate() catch |err| {
        if (err != error.InvalidAtomIndex) return err;
        return;
    };
    return error.ExpectedInvalidAtomIndex;
}
