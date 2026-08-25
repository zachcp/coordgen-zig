# Fuzzing platform and mechanism

Status: decided by `cgz-r15` against the pinned toolchain
`0.17.0-dev.1516+8a4b5424d`. The harness is implemented by `cgz-7v2.4.4`, which
remains open for the gaps named at the end of this document.

The `fuzz` build step is no longer reserved. `zig build fuzz` runs the targets
against their committed corpora and asserts the harness's own preconditions;
`tools/run-fuzz` drives the search, with per-target budgets recorded in
[`conformance/fuzz_budgets.tsv`](../../conformance/fuzz_budgets.tsv).

Targets, and where each lives:

| Surface | Target | Source |
|---|---|---|
| Native | `generation invariants` | `tests/fuzz_targets.zig` |
| Native | `hostile input` | `tests/fuzz_targets.zig` |
| C ABI | `c abi contract` | `src/c_abi/exports.zig` |

The C ABI target sits beside the `export fn` declarations because they are not
`pub`, so no other file can call them without re-declaring the symbols. That
placement is a linkage consequence, not a preference.

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
- A crash prints `error: test '...' terminated with signal ABRT; input saved
  to '.zig-cache/f/crash'` — but **`zig build` still exits 0**. See the
  correction below; this line previously claimed a non-zero exit.
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

## Corrections, 2026-08-24 (cgz-r27)

Claims above were measured against the same pin they were recorded on, and
corrected here rather than edited away, because a harness built on any of them
would have failed silently in the way this document exists to prevent.

One of these corrections is a correction to an earlier correction in this same
section, and that is left visible on purpose: a single probe was generalised
into a rule, and the real harness then contradicted it. Reproduce with
`tools/check-fuzz-seed-promotion --self-test` and with `tools/run-fuzz`.

**Which file holds the reproducing seed is not fixed, and neither source can
be trusted.** Both behaviours occur on this same pin.

In a small standalone probe, three runs from cleared caches:

| source | size | reproduces |
|---|---|---|
| `.zig-cache/f/crash` | 16 bytes, the crashing input, identical every run | yes, 3/3 |
| `.zig-cache/f/<coverage_id>/{0,1}` | 16 bytes, near misses, different every run | no, 0/3 |

That probe is what an earlier revision of this section generalised from, when
it asserted `crash` is always the right source. **That was wrong.** In the real
harness (`tests/fuzz_targets.zig`), with a deliberately planted failure,
`.zig-cache/f/crash` is written **0 bytes** and stays 0 bytes — which is what
`cgz-r27` reported in the first place, and what the small probe simply failed
to reproduce.

Ruled out as the cause, each by measurement:

- **Not input size.** A probe drawing up to 40 values before failing still
  produced a 17-byte `crash`; the fuzzer minimises, so a large draw does not
  mean a large recorded input.
- **Not `-Dfuzz-filter`.** Filtered and unfiltered runs of the real harness
  both produced 0 bytes.

The root cause is not yet identified, and is `cgz-1js`. Until it is, the only
durable rule is the one that does not depend on knowing: **try every candidate
and keep whichever reproduces**, and commit nothing that does not.
`tools/run-fuzz` does exactly that, and refuses to write a seed it could not
replay — it reports the crash and asks for the input to be captured by hand
rather than committing one that tests nothing.
`tools/check-fuzz-seed-promotion` reports which source is winning on a given
toolchain.

**A DISCOVERED crash exits 0, and that is narrower and worse than it sounds.**
Measured for both failure shapes — a returned error and a `@panic` reporting
`signal ABRT`, the exact case the line above described. In each, the crash is
printed, the seed is written, and `zig build test --fuzz=N` exits **0**.

The exit status is not uniformly broken, which is precisely what makes it a
trap. Measured on one target minutes apart:

| plant | when it fires | `zig build --fuzz` exit |
|---|---|---|
| fires on every input | the initial smoke input, before the search | **1** |
| needs discovery | found at run ~3400 | **0** |

So the exit status is correct for failures the fuzzer did not have to find, and
wrong for exactly the ones it exists to find. Anyone spot-checking it with a
trivially-failing target would conclude it works.

This is a third silent failure in the family cgz-r27 names, and the most
dangerous of the three: the backend can be right, the seed can be written, the
crash can be found, and CI is still green. A harness must decide pass/fail from
the run's output and the artifacts it leaves, never from that exit status.

**Why the tool does not hardcode the corrected answer.** `cgz-r24` will move to
released Zig, and either behaviour may change again. The checked-in tool
determines empirically which source yields a reproducing seed and fails when
none does, so the answer is re-measured rather than re-argued.

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

3. **Do not trust `zig build --fuzz`'s exit status, and verify every promoted
   seed by replay.** Both of this item's earlier instructions were measured
   false on the pin and have been corrected below. What survives is the
   requirement they were trying to serve: a promoted seed must be *proven* to
   reproduce before it is committed, and pass/fail must be read from the run's
   output and artifacts rather than from its exit code.

## Known gaps

`cgz-7v2.4.4` stays open for these. Each is a thing the harness does not do,
stated here so that "the fuzz step is green" is not read as more than it is.

1. **Automatic seed promotion does not complete on this toolchain.** The
   reproducing input for a failing target is not reliably persisted anywhere
   the driver can find it (`cgz-1js`). `tools/run-fuzz` detects the failure and
   prints it, tries every candidate, and refuses to commit one that does not
   reproduce - so nothing degrades silently, but a finding must currently be
   captured by hand from the stack trace.

2. **Coverage is not a gate.** The budgets are iteration counts. Nothing
   asserts that a run reached any particular part of the tree, so a target
   that stops exercising a module would still pass its budget. The
   allocation-site ratchet in `conformance/allocation_sites.tsv` does that job
   for the ordinary test suite; the fuzz targets have no equivalent.

3. **Three targets is not "all major input boundaries and algorithm domains",**
   which is what `cgz-7v2.4.4` asks for. There is no target for the residue and
   protein paths, for template loading, or for the `Options` surface.

4. **The C ABI target does not fuzz ownership SEQUENCES.** It exercises one
   generate/free cycle per iteration. Interleaved generate/free orderings
   across several results, which is where a lifetime defect would live, are not
   reached.
