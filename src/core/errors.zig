/// Errors returned by allocating Zig entry points. Assertions are reserved for
/// proven internal invariants; caller-controlled data always returns an error.
pub const Error = error{
    EmptyGraph,
    TooManyItems,
    InvalidAtomicNumber,
    InvalidBondOrder,
    InvalidAtomIndex,
    InvalidStereo,
    InvalidCoordinate,
    InvalidOption,
    InvalidMapping,
    Unsupported,
    OutOfMemory,
};

/// Stable numeric errors crossing the C ABI. Never renumber existing values.
pub const ErrorCode = enum(u32) {
    ok = 0,
    empty_graph = 1,
    too_many_items = 2,
    invalid_atomic_number = 3,
    invalid_bond_order = 4,
    invalid_atom_index = 5,
    invalid_stereo = 6,
    invalid_coordinate = 7,
    invalid_option = 8,
    invalid_mapping = 9,
    out_of_memory = 10,
    unsupported = 11,
    internal = 12,
};

pub fn code(err: Error) ErrorCode {
    return switch (err) {
        error.EmptyGraph => .empty_graph,
        error.TooManyItems => .too_many_items,
        error.InvalidAtomicNumber => .invalid_atomic_number,
        error.InvalidBondOrder => .invalid_bond_order,
        error.InvalidAtomIndex => .invalid_atom_index,
        error.InvalidStereo => .invalid_stereo,
        error.InvalidCoordinate => .invalid_coordinate,
        error.InvalidOption => .invalid_option,
        error.InvalidMapping => .invalid_mapping,
        error.OutOfMemory => .out_of_memory,
        error.Unsupported => .unsupported,
    };
}
