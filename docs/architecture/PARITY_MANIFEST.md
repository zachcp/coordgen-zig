# Corpus stability classification (the parity manifest)

Status: implemented by `cgz-r05` against upstream commit
`d20e735d96480385b2e257522288004038a08cc9`.

This document covers the **stability catalog** half of the parity manifest:
which comparison tier each input and each observable is entitled to, and the
evidence that put it there. The other half — the requirements/coverage matrix
mapping upstream tests, APIs, fixtures, and templates to native tests — is
`cgz-7v2.4.1`. The catalog says *how* to compare; the matrix says *what* must
be covered.

The catalog is consumed under the frozen three-tier policy in
[`COMPARISON_SEMANTICS.md`](COMPARISON_SEMANTICS.md), including its
reflection-safe coordinate normalization and same-build baseline rule.

## Pinned public-test baseline

`tools/check-upstream-test-inventory.py` derives the public-test inventory
directly from the two pinned Boost.Test translation units. It recognizes an
assertion site as a C/C++ token matching
`BOOST_(CHECK|REQUIRE|TEST)[_A-Z]*`, followed by optional whitespace and `(`.
`BOOST_AUTO_TEST_CASE` declarations are test cases, not assertion sites.
Comments and calls made dynamically by a site do not add to the static count.

At `d20e735d96480385b2e257522288004038a08cc9` the reproducible baseline is:

| Item | Count |
|---|---:|
| `test_coordgen.cpp` cases | 16 |
| `test_smilesparser.cpp` cases | 4 |
| Assertion macro sites | **102** |
| Test fixtures | 6 |
| Fixture atoms / bonds | 133 / 142 |
| Example programs | 1 |

The 102 sites are 8 `CHECK`, 1 `CHECK_EQUAL`, 1
`CHECK_EQUAL_COLLECTIONS`, 2 `CHECK_MESSAGE`, 40 `REQUIRE`, 19
`REQUIRE_EQUAL`, and 31 `TEST`. The earlier 95-site audit was not a different
counting convention: it searched for `BOOST_REQUIRE(` literally and missed
seven valid invocations written as `BOOST_REQUIRE (`. The build runs the
token-aware checker as part of `upstream-oracle`, so whitespace cannot silently
reduce the baseline again.

### Corrections to the upstream test system

The pinned harness has six defects. They are evidence to preserve, not
behavior to reproduce:

| Upstream defect | Correct local treatment |
|---|---|
| CTest registers `test_smilesparser` with the `test_coordgen` executable, so its four cases never run. | The inventory reads both source files independently. The requirements matrix and rehosted tests must carry all 20 named cases; CTest discovery is never the source of truth. |
| The example test is registered even when example building is disabled. | The oracle example exists only inside the explicit `-Denable-oracle=true` build branch. Ordinary and oracle-disabled builds neither build nor run it. |
| The example prints coordinates but asserts nothing. | The Zig build requires exit code 0 **and** exact stderr `(-50, 0)  (0, 0)`, turning the observed output into an executable assertion. |
| Test setup fetches maeparser from moving branch `master`. | The production package has no maeparser dependency. The only upstream package is an immutable commit URL plus Zig package hash, cross-checked against `upstream/coordgenlibs.lock`. |
| Valgrind options are configured but no CI job invokes CTest's MemCheck action. | Repo-owned Zig fixture parsing uses `checkAllAllocationFailures`; production tests run without the C++/maeparser harness. Memory claims must come from executed Zig allocation/fuzz gates, never dormant CMake options. |
| `coordgenBasicSMILES.h` uses `sketcherMinimizerAtom` members without including its definition, making compilation depend on include order. | Native corpus construction is Zig and has no textual include order. The conformance-only C++ validator includes each complete atom, bond, and molecule type explicitly before using the upstream helper; that workaround is not exposed as a local header contract. |

The checker deliberately verifies that all six defects are still present in
the pinned source. If a future pin fixes one, the gate fails and requires this
correction table and its corresponding local test to be re-audited instead of
quietly retaining stale assumptions.

## What is measured

The same corpus runs through three oracle builds that differ in exactly one
variable:

| Axis | Baseline | Perturbation |
|---|---|---|
| Architecture | native target | the other CPU architecture of the same OS, cross-compiled |
| Heap address order | platform allocator | `conformance/allocator_order.cpp`: a global `operator new` handing out **descending** addresses |

The allocator axis exists because upstream keys 94 `std::set`/`std::map`
declarations on object pointers, whose iteration order is heap address order.
Inverting that order changes every one of those comparisons without touching a
line of upstream source.

Both perturbations are compared per **member × observable**, never per member
alone. That distinction is the point: coordinates going unstable must not
demote ring sets, fragment trees, or canonical mappings, which is exactly what
`cgz-r13`'s decision record requires.

## Corpus

`src/conformance/corpus.zig` generates both partitions from nothing but a
partition name and an index — no files, no clock, no allocator-dependent
iteration — so every build under comparison sees identical inputs.

- **PRNG**: xorshift64\* (Vigna 2016), pinned by value: seeding constants,
  eight warm-up draws, and modulo reduction *including its bias*. A different
  generator is a different corpus even if it is equally random.
- **`adversarial`**, 2000 members: random valences, free metals (Fe, Zn),
  formal charges, one to three disconnected components, spanning tree plus
  ring closures that produce both ordinary rings and macrocycles. This is the
  population behind the divergence numbers below.
- **`drug_like`**, 7 members: two fused polycyclics, a steroid skeleton, a
  seventeen-membered macrocycle, a bridged bicyclic, a peptide-like chain, and
  a four-ring biaryl chain. Committed as tables, and checked at build time
  against upstream's own SMILES parser (`conformance/smiles_reference.cpp`) by
  diffing two independently produced dumps.

Both the corpus size and the partition list are part of the contract: a
classification is only meaningful for the population it was measured on.

## Results

Measured on aarch64-macos against x86_64-macos, Zig
`0.17.0-dev.1516+8a4b5424d`, oracle built `ReleaseFast`. Counts are members
whose observable differed under that axis.

### Stable under both perturbations

Zero divergent members out of 2000 adversarial and 7 drug-like, on both axes:

`status`, `clean_pose`, `probe_status`, `probe_clean_pose`,
`input_to_internal`, `internal_to_input`, `morgan_ranks`,
`effective_bond_orders`, `bond_displays`, `atom_stereo`, `rings`,
`rings_set`, `fragments`, `fragments_set`, `template_mappings`,
`template_mappings_set`, `components`, `components_set`.

### Divergent

| Partition | Observable | Members | Architecture | Allocator order |
|---|---|---:|---:|---:|
| adversarial | `coordinates` | 2000 | 1414 (70.7%) | 1 (0.05%) |
| adversarial | `component_transforms` | 2000 | 1972 (98.6%) | 4 (0.2%) |
| adversarial | `dofs` | 2000 | 36 (1.8%) | 0 |
| adversarial | `dofs_set` | 2000 | 33 (1.7%) | 0 |
| adversarial | `dof_penalties` | 2000 | 19 (1.0%) | 0 |
| drug_like | `component_transforms` | 7 | 7 (100%) | 0 |
| drug_like | everything else | 7 | 0 | 0 |

Largest coordinate deviation over the adversarial partition: **49.4 bond
lengths** across architectures, **0.58 bond lengths** across allocator order.
Size buckets present: 187 tiny (<10 atoms), 764 small, 674 medium, 332 large,
50 huge (≥70 atoms).

The 0.58 figure supersedes the *739 units, roughly 15 bond lengths* recorded in
`cgz-r13` and in the epic's summary of it, which came from the standalone C++
review probe. This pipeline measures the same single member 25× smaller. The
correction matters because that number is the whole empirical case for the
parity ceiling: the one input a typed-index port structurally cannot match is
off by half a bond length, not by a relayout.

The adversarial coordinate counts — 1414/2000 across architectures and 1/2000
across allocator order — reproduce the review-probe measurements recorded in
`cgz-r05` through an entirely different pipeline: a Zig corpus generator, the
C ABI oracle adapter, and this classifier, rather than a standalone C++ probe.

`component_transforms` diverges far more often than the coordinates it is
derived from — including on all seven drug-like members, whose coordinates are
bit-identical. Inspecting the payloads shows differences of one to two units
in the last place of the transform matrix, not different placements. It is a
derived float observable and belongs to a tolerant tier, never an exact one.

## How to read this

- **Structural observables are stable under both perturbations.** Ring sets
  and order, fragment trees and flags, canonical and input mappings, morgan
  ranks, template mappings, effective bond orders, and the `clean_pose` flag
  did not move on any of the 2000 adversarial members under either axis.
  Those observables are entitled to exact comparison even for inputs whose
  coordinates are not.
- **Coordinates are the unstable observable.** They are the reason the epic's
  two-tier "exact or tolerant" split is insufficient.
- **Derived floats are reported separately from the structures that carry
  them.** DOF penalties and component transforms are emitted as their own
  observables precisely because a one-unit-in-the-last-place difference in a
  derived float is a different finding from a DOF that changed state or a
  component that gained an atom.

## What is gated, and what is published

These are deliberately different artifacts.

- `conformance/parity_expectations.tsv` is the **gate**. It records, per
  partition and observable, the fraction of members allowed to fail each axis.
  An observable with no recorded expectation fails the gate: a new observable
  must be classified before anything can depend on it. This file is portable
  because it states claims that hold on any host.
- `conformance/parity_manifest.tsv` is the **published evidence**: the
  per-member classification and the coordinate deviations that produced the
  numbers above. Coordinate values and deviations are
  per-(architecture, toolchain, optimize-mode) artifacts, so this file
  describes the build that produced it and is never a cross-platform fixture.
  Regenerate it with `zig build parity-manifest -Denable-oracle=true`.

Running the gate:

```sh
zig build corpus-check -Denable-oracle=true
```

The architecture axis needs a host that can execute the other architecture's
binaries — Rosetta on aarch64 macOS, `-fqemu` with `qemu-user` on Linux. If it
cannot, `corpus-check` fails with that message rather than skipping the axis:
a skipped Run step still leaves a build green, and an axis that quietly
disappears from a classification is worse than a red build.

## Notes for the beads that consume this

- `cgz-r06` (comparison semantics): the per-observable columns are the input
  to tier assignment. Structural observables stable on both axes qualify for
  the exact tier on the whole corpus, not only on the stable partition.
- `cgz-r13` (parity ceiling): settled in
  [COMPATIBILITY_POLICY.md](COMPATIBILITY_POLICY.md). The ceiling is per
  observable, is scoped to the **allocator-order axis only**, and is the
  enumerated `order_unstable` column of `parity_manifest.tsv` rather than any
  fraction in `parity_expectations.tsv`. At the pin it is four member ×
  observable groups: `adversarial/917`, `1538`, `1588` (component transforms
  only, coordinates bit-identical) and `adversarial/1695` (coordinates, 0.580
  bond lengths). The 1414/2000 architecture figure is **not** a ceiling — the
  same-build rule means both sides of a differential comparison share an
  architecture.
- `cgz-7v2.4.2` (differential runner): reuse `tests/oracle_corpus_run.zig`'s
  dump format. Oracle and native must run in the same process, same build,
  same target.
