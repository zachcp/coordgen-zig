pub const api = @import("api");
pub const Input = api.Input;
pub const AtomInput = api.AtomInput;
pub const BondInput = api.BondInput;
pub const ResidueInput = api.ResidueInput;
pub const ResidueInteractionInput = api.ResidueInteractionInput;
pub const Options = api.Options;
pub const Result = api.Result;
pub const bond_length = api.bond_length;

// `test {}` blocks are only analyzed by `zig build test`, never by a plain
// `zig build`/install, and merely naming a module in a `pub const` does not
// force the compiler to walk its declarations either - verified empirically:
// a `pub const c_abi = @import("c_abi");` with nothing else referencing it
// produced an installed libcoordgen.a with zero symbols from that file. Only
// a `comptime` block that actually references the import forces the
// compiler to discover its `export fn`s, so that reference has to live on
// *this* real top-level path, not merely inside the test block below, or
// coordgen_generate/coordgen_result_free would silently be absent from the
// installed archive while still compiling and passing in test builds.
pub const c_abi = @import("c_abi");
comptime {
    _ = c_abi;
}

test {
    // Internal contract tests are aggregated only in test builds. Production
    // consumers see the safe public API above, not mutable model internals or
    // conformance types.
    _ = api;
    _ = c_abi;
    _ = @import("core");
    _ = @import("model");
    _ = @import("geometry");
    _ = @import("module_layers");
}
