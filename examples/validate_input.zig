const coordgen = @import("coordgen");

pub fn main() !void {
    const atoms = [_]coordgen.AtomInput{
        .{},
        .{ .atomic_number = .oxygen },
    };
    const bonds = [_]coordgen.BondInput{
        .{ .start = 0, .end = 1, .order = .double },
    };
    const input: coordgen.Input = .{ .atoms = &atoms, .bonds = &bonds };
    try input.validate();

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
