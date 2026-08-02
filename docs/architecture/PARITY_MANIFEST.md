# Corpus stability classification (the parity manifest)

Status: implemented by `cgz-r05` against upstream commit
`d20e735d96480385b2e257522288004038a08cc9`.

This document covers the **stability catalog** half of the parity manifest:
which comparison tier each input and each observable is entitled to, and the
evidence that put it there. The other half — the requirements/coverage matrix
mapping upstream tests, APIs, fixtures, and templates to native tests — is
`cgz-7v2.4.1`. The catalog says *how* to compare; the matrix says *what* must
be covered.

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
- `cgz-r13` (parity ceiling): the ceiling is per observable. Coordinates on
  order-unstable members are the only observable this corpus shows to be
  unreproducible in the oracle itself.
- `cgz-7v2.4.2` (differential runner): reuse `tests/oracle_corpus_run.zig`'s
  dump format. Oracle and native must run in the same process, same build,
  same target.
