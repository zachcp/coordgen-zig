# Fuzzing platform and mechanism

Status: decided by `cgz-r15` against the pinned toolchain
`0.17.0-dev.1516+8a4b5424d`. Implemented by `cgz-7v2.4.4`; the iteration
budgets and pass conditions are in [`SUCCESS_CRITERIA.md`](SUCCESS_CRITERIA.md).

The `fuzz` build step is reserved and fails with its owning bead until that
implementation lands, per the `cgz-r21` no-relaxation policy.

## Decision

Fuzz **both CI platforms** — `ubuntu-24.04` (x86\_64, ELF) and `macos-14`
(aarch64, Mach-O) — with **one engine**, Zig's built-in coverage-guided fuzzer,
over **two layered surfaces**.

| Surface | Targets | Why it needs its own targets |
|---|---|---|
| Native modules | parsing, graph validation, geometry, topology, macrocycle state, bounded full generation | `Smith` builds valid-but-hostile deep algorithmic states that random bytes never reach |
| Stable C ABI (`include/coordgen_abi.h`) | malformed input, ownership and lifetime sequences, error paths | this is the shipped surface; its ownership rules are not reachable from native module targets |

Both surfaces use `std.testing.fuzz` with `Smith`. The C ABI surface is a Zig
fuzz target that *calls* the exported C entry points as a `Smith`-driven state
machine; it is not an external AFL++/libFuzzer harness.

## Why not AFL++ / libFuzzer

The pinned compiler does not support libFuzzer instrumentation:

```
$ zig build-exe -fsanitize=fuzzer probe.zig
error: unrecognized parameter: -fsanitize=fuzzer
```

An AFL++/libFuzzer harness would therefore require a **second, unchecksummed C
toolchain installed and version-managed on both runners**. `cgz-r01` exists to
make the build reproducible from one pinned, checksummed toolchain, and that
cost would only be justified if the built-in fuzzer did not work. It does work,
on both platforms, and was measured doing so.

Using one engine also keeps a single seed format, so a failure found on either
surface is promoted and replayed by the same mechanism.

## Toolchain status, re-checked

`zig build --fuzz` was broken on macOS aarch64 in Zig 0.16.0: the shipped test
runner failed to compile in fuzz mode at `lib/compiler/test_runner.zig:566`
(`*const debug.StackTrace` versus `*builtin.StackTrace`). **That blocker is gone
on the pin** — the file no longer references `StackTrace`, and fuzz-mode builds
succeed.

Measured on the pin, host aarch64-macos, each run from a cleared `.zig-cache`:

- A five-deep coverage staircase — a `@panic` reachable only after five separate
  byte comparisons — was found in **75 and 78 runs** on two independent runs.
  Unique runs climbed 0 → 13 and 0 → 15. That is coverage-guided progress, not
  random search.
- A crash **fails the build**: non-zero exit,
  `error: test '...' terminated with signal ABRT`.
- A fuzz-instrumented test binary **cross-compiles cleanly** for
  `x86_64-linux-gnu`.
- Throughput anchor: 1,000,007 iterations of a trivial target in **15.1 s** wall
  clock, single instance.

`lib/fuzzer.zig` supports exactly two object formats, `.elf` and `.macho`; any
other is a `@compileError`. Both CI platforms are covered. A Windows/COFF leg
could never run this gate, so adding one would require revisiting this decision.

## Regression-seed promotion

Every minimized failure is committed as a byte seed under version control and
replayed through `FuzzInputOptions.corpus`:

```zig
test "…" {
    try std.testing.fuzz({}, testOne, .{ .corpus = promoted_seeds });
}
```

In an ordinary `zig build test` — no fuzzer engine, no `--fuzz` — `std.testing.fuzz`
feeds each corpus entry through a `Smith` backed by those bytes, plus an empty
input as a smoke test. Regression coverage therefore runs **on every platform
and every build**, including builds and platforms where the engine never runs.

This was verified by falsification rather than assumed: a probe that returns an
error on any non-empty decode does fail under plain `zig build test`, so the
seeds are genuinely executed and not silently ignored.

Promoted seeds are byte strings interpreted by a specific target's `Smith` call
sequence. They are **coupled to that target's decode order**: changing the order
or type of `Smith` draws in a target reinterprets its existing seeds. A target's
seeds must be regenerated when its decode sequence changes.

## Operating constraints

These are properties of the pinned toolchain, each measured. The harness must
respect them.

1. **Always pass an explicit iteration limit.** Bare `--fuzz` implies `--webui`
   and blocks forever; `--fuzz=<N>` (with optional `K`/`M`/`G` suffix) runs
   headless and bounded. Iteration limits, not wall-clock, are the budget unit —
   they are reproducible across runners of differing speed.

2. **Assert the compiler backend.** `test_runner.zig` sets `need_simple` for the
   `stage2_aarch64`, `stage2_powerpc`, and `stage2_riscv64` backends, and
   `std.testing.fuzz` then **returns immediately**. The pin currently selects
   `stage2_llvm` on aarch64-macos, so fuzzing is live — but if a self-hosted
   backend becomes the default, every fuzz target silently becomes a no-op that
   still reports green. That is precisely the `cgz-r21` false-green failure
   mode, so the harness must assert `builtin.zig_backend` is not a `need_simple`
   backend instead of trusting the default.

3. **Do not promote seeds from the fuzzer's `crash` file.** On the pin,
   `.zig-cache/f/crash` is written **0 bytes** even when the crashing input is
   known — reproduced on two independent crashing runs from cleared caches. The
   persistent corpus directory `.zig-cache/f/<coverage_id>/{0,1,2,…}` does hold
   real, replayable seeds; promotion must read that directory.
