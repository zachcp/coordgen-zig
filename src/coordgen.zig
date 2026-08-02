pub const api = @import("api");
pub const Input = api.Input;
pub const AtomInput = api.AtomInput;
pub const BondInput = api.BondInput;
pub const ResidueInput = api.ResidueInput;
pub const ResidueInteractionInput = api.ResidueInteractionInput;
pub const Options = api.Options;
pub const Result = api.Result;
pub const bond_length = api.bond_length;

test {
    // Internal contract tests are aggregated only in test builds. Production
    // consumers see the safe public API above, not mutable model internals or
    // conformance types.
    _ = api;
    _ = @import("core");
    _ = @import("model");
    _ = @import("geometry");
    _ = @import("c_abi");
    _ = @import("module_layers");
}
