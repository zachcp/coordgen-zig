const std = @import("std");

// Smith's deterministic replay consumes one little-endian u64 per value. This
// describes a validator-clean ten-carbon chain followed by four generates and
// four frees across distinct result slots. Both supported CI architectures are
// little-endian; make that assumption explicit rather than hiding it in an
// opaque escaped byte string.
comptime {
    if (@import("builtin").target.cpu.arch.endian() != .little) {
        @compileError("C ABI bootstrap corpus requires little-endian Smith words");
    }
}

const atom = [_]u64{ 6, 128, 0 }; // carbon, zero charge, unspecified stereo
const bond = [_]u64{1}; // single
const operations = [_]u64{
    8, // operation count
    0, 1, 1, 1, 2, 1, 3, 1, // generate slots 0..3
    0, 0, 1, 0, 2, 0, 3, 0, // free slots 0..3
};
const words = [_]u64{ 1, 10 } ++
    atom ++ atom ++ atom ++ atom ++ atom ++ atom ++ atom ++ atom ++ atom ++ atom ++
    bond ++ bond ++ bond ++ bond ++ bond ++ bond ++ bond ++ bond ++ bond ++
    operations;

/// Deterministic Smith input that starts from a valid generated molecule and
/// reaches a multi-result ownership sequence from an empty fuzzer cache.
pub const valid_sequence: []const u8 = std.mem.asBytes(&words);
