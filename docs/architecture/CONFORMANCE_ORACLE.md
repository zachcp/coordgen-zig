# Conformance oracle and fixture reading

Status: implemented by `cgz-r03` against upstream commit
`d20e735d96480385b2e257522288004038a08cc9`. Scope is the oracle *build* and
fixture *reading*. Differential comparison against the native port is
`cgz-7v2.4.2`; the oracle-facing ABI split is `cgz-r04`.

## The oracle is built by Zig, not by CMake

`zig build upstream-oracle -Denable-oracle=true` compiles the fourteen pinned
upstream translation units into a static library through the Zig build graph:

- `-std=c++17 -DIN_COORDGEN`, `link_libcpp`, no CMake, no system toolchain.
- The file list in `build.zig` is written out explicitly. An upstream file
  appearing or disappearing is then a build failure, not a silent change in
  what the oracle contains.
- The library is **never installed**. `tools/check-install-isolation` fails if
  an oracle-named file or an upstream C++ symbol reaches an install prefix.
- The dependency is lazy (`.lazy = true` plus `b.dependencyLazy`) and gated
  behind `-Denable-oracle=true`, so a build that never asks for conformance
  never fetches upstream sources. CI proves this by running `zig build test`
  with an empty cache and no network.

Upstream's example program (`example_dir/example.cpp`) is compiled against the
static library and run as part of the step; its stderr must be exactly
`(-50, 0)  (0, 0)` — its two atoms one bond length (50) apart on the x axis.
Upstream's README points at the example but does not publish its output, so
that string is a baseline recorded here against the pin, not an upstream
promise. It is the cheapest end-to-end proof that the library links and runs.

### Do not try to build upstream's CMake test target

The upstream **library** is Boost-free and maeparser-free. The upstream
**test target** is not: it needs Boost.Test, Boost.filesystem, and maeparser.
Building it through Zig is not a matter of more flags, and attempting it is a
week lost. Re-hosting the harness and the fixture reading on the Zig side is
what makes the oracle buildable at all — which is why the minimal MAE reader
below is a prerequisite, not a convenience.

## Minimal MAE reader

`src/conformance/mae.zig` reads exactly what upstream's
`sketcherMaeReading.h` consumes, and nothing else:

| Block | Fields |
|---|---|
| `m_atom[N]` | `i_m_atomic_number`, `r_m_x_coord`, `r_m_y_coord` |
| `m_bond[N]` | `i_m_from`, `i_m_to`, `i_m_order` |

Contract:

- Fields are resolved **by name**, never by column. The pinned fixtures put
  `i_m_atomic_number` at four different column positions, and one
  `templates.mae` variant inserts `i_m_formal_charge` ahead of the last
  column.
- Maestro atom references are 1-based; `Bond.from`/`Bond.to` are 0-based
  positions into `Structure.atoms`, converted once, in the reader.
- Coordinates are kept as `f64`. Upstream narrows to `float` when it builds
  `sketcherMinimizerPointF`; that narrowing belongs to the oracle adapter, not
  to the file reader.
- Every other block (`f_m_ct` scalar properties, `m_depend[N]`, the
  `s_m_m2io_version` header) is checked structurally and discarded.
- One `Reading` owns all its arrays through one arena and requires exactly one
  `deinit`. Allocation failure is reported, never silently truncated; the
  reader is exercised under `checkAllAllocationFailures`.

Syntax the pinned fixtures actually exercise, all covered by tests:

- A quoted value in `test/macrocycle.mae` spans twenty-one physical lines, and
  its continuation lines start with digits at column 0. A line-oriented reader
  desynchronizes by twenty lines here and then mistakes `m_atom` for a value.
- `f_m_ct` property names contain backslashes and parentheses, for example
  `s_sd_PUBCHEM\_CACTVS\_SUBSKEYS_(undefined)`.
- Data rows carry the `<>` missing-value sentinel in columns this reader does
  not consume, and quoted values (`""`, `"N6"`) inside `templates.mae` rows.
- `# … #` comments, runs of two or more spaces between columns, trailing
  spaces on block delimiters, and `-0.000000` coordinates.

Deliberate deviations from upstream behaviour, all of them stricter:

| Situation | Upstream | This reader |
|---|---|---|
| Bond referencing a nonexistent atom | `catch (std::out_of_range)` around the whole bond loop silently drops the remaining bonds | `error.AtomIndexOutOfRange` |
| Row index column out of sequence | ignored | `error.InvalidRowIndex` |
| Missing required field | maeparser throws | `error.MissingRequiredProperty` |

Parse failures carry a `Diagnostic` with the line number and offending token.

## Fixture expectations

`build.zig` holds a `fixture_expectations` table covering all six upstream
test fixtures plus `templates.mae`, and `tests/mae_fixture_check.zig` reads
each file from the lazy dependency and compares totals.

| Fixture | Structures | Atoms | Bonds |
|---|---:|---:|---:|
| `test/test.mae` | 1 | 26 | 26 |
| `test/test_mol.mae` | 1 | 27 | 31 |
| `test/testChirality.mae` | 1 | 9 | 8 |
| `test/macrocycle.mae` | 1 | 59 | 67 |
| `test/metalZobs.mae` | 1 | 6 | 5 |
| `test/nonterminalMetalZobs.mae` | 1 | 6 | 5 |
| `templates.mae` | 82 | 1704 | 1963 |

`metalZobs.mae` and `nonterminalMetalZobs.mae` have identical totals by
design: they differ only in `f_m_ct` provenance properties
(`s_m_Source_Path`, `s_m_Source_File`, `i_m_Source_File_Index`), not in any of
the six fields.

The totals — including atomic-number, bond-order, and coordinate sums — are
derived by `tools/mae-fixture-totals.py`, a second implementation that
tokenizes the fixtures itself. The expectations are therefore not a digest of
this repository's own output: a mismatch means the two readers disagree. To
re-derive them after a pin change:

```sh
zig build --fetch=all   # populates the lazy oracle package
tools/mae-fixture-totals.py "$ORACLE"/test/*.mae "$ORACLE"/templates.mae
```

`$ORACLE` is the fetched package directory. Copy the emitted `KEY=VALUE` lists
into `fixture_expectations`.

## What the step does and does not prove today

Proven: the pinned sources compile and link through Zig on `aarch64-macos` and
cross-compile unchanged to `x86_64-macos`; the resulting library reproduces
upstream's documented example output; all seven `.mae` files parse to totals
that an independent implementation agrees with.

Since then `cgz-r04` added the conformance probe ABI over the same oracle, and
`cgz-r05` added the corpus stability classification — see
[PARITY_MANIFEST.md](PARITY_MANIFEST.md). `zig build conformance` now runs the
oracle step plus `corpus-check`.

Not yet proven, and owned elsewhere: any comparison of oracle output against
the native port (`cgz-7v2.4.2`).
