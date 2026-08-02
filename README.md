# coordgen-zig

Native Zig port of CoordGen, currently at the conformance-foundation stage.
The production `coordgen` module is standard-library-only. The pinned C++
implementation is an opt-in test oracle, never a production dependency.

Bootstrap the exact compiler before invoking Zig:

```sh
ZIG=$(tools/bootstrap-zig)
"$ZIG" build test
"$ZIG" build package-check
```

Fetch and validate the conformance dependency explicitly:

```sh
"$ZIG" build --fetch=all
"$ZIG" build -Denable-oracle=true upstream-oracle
tools/verify-upstream
```

See `docs/TOOLCHAIN.md`, `docs/PROVENANCE.md`, and `docs/DEPENDENCIES.md`
before changing compiler or dependency pins.
