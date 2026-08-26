# Quantified success criteria

Status: set by `cgz-r14`. Every criterion below is either a number backed by a
measurement in this repo, or an explicitly named gap with a trigger that fixes
it. Nothing here may be relaxed to make a gate pass (`cgz-r21`).

Consumers: `cgz-7v2.4.3` (property/malformed/OOM/concurrency),
`cgz-7v2.4.4` (fuzzing), `cgz-7v2.4.7` (performance and release gates).

## Two different comparisons, never conflated

This is the most common way to misread the numbers below.

| | Compares | Gated by | What it bounds |
|---|---|---|---|
| **Stability classification** | oracle vs *perturbed oracle* (other arch, reversed heap order) | `conformance/parity_expectations.tsv` | how much the **upstream algorithm** moves under perturbation |
| **Parity gate** | **native** vs oracle, same process/build/target | this document | whether **the port** reproduces upstream |

`parity_expectations.tsv` is *not* a native tolerance. It says which tier each
input × observable is entitled to. The tolerances for native-vs-oracle are here.

## Corpus

Generated deterministically by `src/conformance/corpus.zig` from a partition
name and index — no files, no clock, no allocator-dependent iteration.

| Partition | Members | `ok` (oracle) | `clean_pose` true | Size buckets |
|---|---:|---|---:|---|
| `adversarial` | 2000 | 2000/2000 (100%) | 348 (17.4%) | 186 tiny, 758 small, 674 medium, 332 large, 50 huge |
| `drug_like` | 7 | 7/7 (100%) | 6 (85.7%) | 1 tiny, 6 small |

Corpus size and partition list are part of the contract. Narrowing either is a
gate relaxation and requires a decision bead.

## Tier 1 — exact

These 18 observables showed **zero divergent members out of 2007** on both
perturbation axes, so they are entitled to exact comparison on the whole corpus:

`status`, `clean_pose`, `probe_status`, `probe_clean_pose`, `input_to_internal`,
`internal_to_input`, `morgan_ranks`, `effective_bond_orders`, `bond_displays`,
`atom_stereo`, `rings`, `rings_set`, `fragments`, `fragments_set`,
`template_mappings`, `template_mappings_set`, `components`, `components_set`.

**Pass condition: 2007/2007 members, byte/value equality, zero mismatches.**
Not a rate — any single mismatch fails. `ok` and `clean_pose` are included:
the oracle achieves `ok` on 100% of both partitions, so native must too.

## Tier 2 — tolerant

Applies to coordinates and derived floats on input × observable pairs the
catalog marks stable. Compared after translation/rotation normalization;
reflection is rejected unless the fixture carries an achirality proof.

**Initial tolerance T = 0.1 bond lengths (5.0 CoordGen units).**

Derived from the measured architecture-axis deviation distribution over the 1414
divergent adversarial members, which is strongly bimodal:

| Deviation (bond lengths) | Members | |
|---|---:|---|
| < 0.05 | 1024 | 72.4% — float jitter, same layout |
| [0.05, 0.1) | 23 | |
| **[0.1, 1.0)** | **33** | **2.3% — sparse valley** |
| [1.0, 2.0) | 57 | |
| ≥ 2.0 | 277 | 23.6% total ≥ 1.0 — a different local minimum |

Only 2.3% of divergent members land between 0.1 and 1.0 bond lengths. That
valley separates "same layout, floating-point jitter" from "the discrete search
found a different minimum", and 0.1 is its strict end.

**This is a prior, not a native measurement.** It bounds the *upstream
algorithm's own* float sensitivity, because no native implementation exists yet
and no native-vs-oracle deviation has ever been measured. T must be
**recalibrated against the first native-vs-oracle baseline** (`cgz-7v2.4.2`).
Recalibration may lower T. Raising it above 0.1 requires a decision bead
recording what is being surrendered, because above this valley a "tolerant"
pass no longer means the layout matched.

## Tier 3 — invariant / statistical

For input × observable pairs the oracle cannot reproduce under perturbation.
Aggregated per partition through `comparison.noWorseThan`.

**This tier is unpopulated at the pin** (`cgz-r30`, applying `cgz-r13`).

An earlier draft of this section read "on this corpus, coordinates on the 1414
arch-divergent adversarial members". That is withdrawn. It was written the same
afternoon as `cgz-r13`, on a parallel branch, and merged fifteen seconds after
it; the two were never reconciled, so it is a merge race rather than a
considered disagreement.

Two independent reasons it cannot stand:

1. **The runner cannot observe that axis.** Oracle and native execute in one
   process, build, and target, enforced by `RunIdentity.requireSameBuild()`.
   Both sides of every differential comparison therefore share an architecture,
   so architecture divergence is not a difference the runner can see, let alone
   excuse. This is structural, not a policy preference.
2. **It would surrender the observable the port exists to produce.** 1414 of
   2000 is 70.7% of the adversarial corpus. Demoting that to "no worse than the
   oracle on aggregate quality metrics" would mean the port never proves its
   coordinates match on the majority of its own corpus.

The ceiling that *is* real is the heap-address-order axis, enumerated per
(member, observable) in [`conformance/parity_ceiling.tsv`](../../conformance/parity_ceiling.tsv)
(`cgz-r26`) — nine rows at the pin. Those members are held to a **tolerant**
comparison against a published per-member bound, not excused to statistics:
`adversarial/1695`'s 0.580 bond lengths sits inside a tolerant band, which is
precisely why `cgz-r13` left this tier empty.

The margins below remain the mechanism if a future measurement ever populates
this tier. Nothing populates it today.

**All five margins in `QualityMargins` stay at their zero defaults**:
`clash_score`, `bond_length_rms`, `bond_angle_deviation`, `bond_crossings`,
`atoms_inside_rings`. Native must be no worse than the same-build oracle on
every metric.

Zero is the correct starting value, not a placeholder: it is the strictest
setting the mechanism allows, and the values currently in `comparison.zig` are
test fixture data, not decided margins. Each nonzero margin requires its own
decision bead naming the metric, the value, and the evidence — so a margin can
never be widened silently to absorb a regression.

## Phase 3 gate, restated measurably

Replaces "ordinary small molecules should produce deterministic usable
coordinates". On the **`drug_like` partition (7 members)**:

1. Native and oracle agree **exactly** on all 18 Tier-1 observables — ring sets,
   canonical ordering, fragment trees, template matches, `ok`, `clean_pose`:
   **7/7, zero mismatches**.
2. Coordinates agree within **T** after normalization: **7/7**. The measurement
   supports the strict form here — all 7 drug-like members are coordinate-stable
   on both axes with max deviation 0.000 bond lengths, so this partition has no
   invariant-tier escape hatch.
3. `ok` on **7/7**, matching the oracle's measured 100%.

`component_transforms` is excluded from (1): it is arch-divergent on 7/7 by
measurement, differing by one to two units in the last place. It is a derived
float and is compared in Tier 2.

## Determinism gates

Five gates, each a hard pass/fail. All run on both CI platforms.

| Gate | Procedure | Pass condition | Catches |
|---|---|---|---|
| Same input twice | run identical input twice in one process | byte-identical results, all observables | static/global state (`cgz-r11`, `cgz-r12`) |
| Two threads | two contexts generating concurrently | byte-identical to serial results; clean under TSan | template race (`cgz-r12`) |
| Reversed atom order | submit input with atom and bond order reversed | Tier-1 observables identical after mapping through `input_to_internal`/`internal_to_input`; coordinates within T | storage order leaking into output |
| Allocator order | descending-address allocator (`conformance/allocator_order.cpp`) | within `parity_expectations.tsv` order ceilings | pointer-keyed container iteration |
| Cross-target | other CPU architecture, cross-compiled | within `parity_expectations.tsv` arch ceilings | float evaluation order |

Reversed atom order is the cheapest and strongest of the five: canonical
ordering is supposed to make it exactly invariant, so it needs no tolerance on
Tier-1 observables at all.

The last two run against the **recorded ceilings**, not against zero, because
the oracle itself is not stable under them. The first three admit no tolerance
on Tier-1 observables.

## OOM and leak gates

- **Leak detection**: `std.testing.allocator` on **every** test. Any leak fails.
- **OOM injection**: `std.testing.FailingAllocator` over **every reachable
  allocation index** for bounded representative inputs to **every allocating
  public entry point** — Zig `generate` and each allocating `coordgen_*` C ABI
  function. Every injected failure must produce a clean error return with no
  leak and no partial-result escape.
- **Coverage counter**: the suite records the count of **distinct allocation
  sites exercised per public entry point**, committed in-repo. The gate fails if
  any entry point's count **decreases**.

The counter exists because "representative inputs" is otherwise unfalsifiable:
without it, the criterion silently weakens as the code grows and new allocation
sites are simply never reached. A decrease is a real regression in coverage,
even when every test still passes.

## Fuzz gates

Platform and mechanism are decided in `cgz-r15`; see
[`FUZZING.md`](FUZZING.md). Budgets:

- **Both CI platforms**, `ubuntu-24.04` and `macos-14`. Neither is exempt.
- **Budgets are iteration limits** (`--fuzz=<N>`), not wall-clock: they are
  reproducible across runners of differing speed. Bare `--fuzz` implies
  `--webui` and hangs, so a limit is mandatory.
- **Per-target limits are recorded in-repo**, so lowering one is a visible diff
  rather than an invisible edit (`cgz-r21`).
- **Ceiling: the whole fuzz step completes within 10 minutes per platform.**
  Per-target limits are calibrated to fit that ceiling on the slower runner.
- **Any crash, panic, leak, or error-log output fails the build.**
- **Every minimized failure is promoted** to a committed byte seed and replayed
  through `FuzzInputOptions.corpus` in the ordinary `zig build test`. Seed count
  only ever grows; removing a seed requires a decision bead.
- **The harness asserts `builtin.zig_backend` is not a `need_simple` backend**,
  because `std.testing.fuzz` returns immediately on those and every target would
  report green while testing nothing.

The one throughput datum available is an anchor, not a budget: 1,000,007
iterations of a *trivial* target in 15.1 s on aarch64-macos. Coordinate
generation is orders of magnitude heavier per iteration, so per-target limits
cannot be extrapolated from it and must be measured when the harness exists.

## Performance gate — thresholds set from the first native baseline

`zig build performance-baseline -Denable-oracle=true` times **native and
oracle in one process**, on the same members, in the same build, target and
optimize mode. `zig build performance-check -Denable-oracle=true` runs the same
binary with `--enforce`, so the gate reads exactly the numbers the baseline
printed. The 100-member population is in `tests/oracle_benchmark.zig`; changing
it requires a decision Bead.

The trigger named here — set the threshold when the first representative
native-vs-oracle baseline exists — has been met, and the thresholds are in
[`conformance/performance_thresholds.tsv`](../../conformance/performance_thresholds.tsv)
with their provenance.

**Method**, unchanged from the frozen decision, plus one asymmetry now worth
stating: the oracle is timed through the C entry points it provides and native
through its Zig API, because the linked oracle already defines the `coordgen_*`
symbols and a second definition is the `cgz-r28` duplicate-symbol defect. The
difference is native's DTO conversion, charged to the oracle side. The ratios
are therefore mildly generous to native.

**A member native cannot lay out is skipped on BOTH sides** and counted. A
ratio between medians taken over different member subsets is not a
like-for-like ratio, and dropping native's hard members while keeping the
oracle's would flatter native by exactly the amount that matters. Coverage is
gated separately, by `min_compared_members`, so losing domain is a failure even
though it leaves no ratio to exceed.

**What the first baseline actually showed**, aarch64-macos, ReleaseFast:

| bucket | median ratio | p95 ratio | members compared |
|---|---|---|---|
| tiny | 0.900–1.014 | 10.2–13.9 | 19 of 20 |
| small | 2.996–3.071 | 3.457–3.483 | 17 of 20 |
| medium | — | — | **0 of 20** |
| large | — | — | **0 of 20** |
| huge | — | — | **0 of 20** |

Two facts in that table matter more than the thresholds:

- **Native cannot lay out any member of the three largest buckets.** The
  performance gate covers 36 of the benchmark's 100 members. That is a
  statement about the port's domain coverage, not about its speed, and it is
  ratcheted so it cannot quietly get smaller.
- **Native is at parity on the tiny median and an order of magnitude slower at
  the tiny p95**, reproducibly. That is `cgz-twy`. A median-only comparison
  would have hidden it, which is why this document requires p95 per bucket.

## OOM and leak gates

- **Leak detection**: `std.testing.allocator` on **every** test. Any leak fails.
- **OOM injection**: `std.testing.FailingAllocator` over **every reachable
  allocation index** for bounded representative inputs to **every allocating
  public entry point** — Zig `generate` and each allocating `coordgen_*` C ABI
  function. Every injected failure must produce a clean error return with no
  leak and no partial-result escape.
- **Coverage counter**: the suite records the count of **distinct allocation
  sites exercised per public entry point**, committed in-repo. The gate fails if
  any entry point's count **decreases**.

The counter exists because "representative inputs" is otherwise unfalsifiable:
without it, the criterion silently weakens as the code grows and new allocation
sites are simply never reached. A decrease is a real regression in coverage,
even when every test still passes.

## Fuzz gates

Platform and mechanism are decided in `cgz-r15`; see
[`FUZZING.md`](FUZZING.md). Budgets:

- **Both CI platforms**, `ubuntu-24.04` and `macos-14`. Neither is exempt.
- **Budgets are iteration limits** (`--fuzz=<N>`), not wall-clock: they are
  reproducible across runners of differing speed. Bare `--fuzz` implies
  `--webui` and hangs, so a limit is mandatory.
- **Per-target limits are recorded in-repo**, so lowering one is a visible diff
  rather than an invisible edit (`cgz-r21`).
- **Ceiling: the whole fuzz step completes within 10 minutes per platform.**
  Per-target limits are calibrated to fit that ceiling on the slower runner.
- **Any crash, panic, leak, or error-log output fails the build.**
- **Every minimized failure is promoted** to a committed byte seed and replayed
  through `FuzzInputOptions.corpus` in the ordinary `zig build test`. Seed count
  only ever grows; removing a seed requires a decision bead.
- **The harness asserts `builtin.zig_backend` is not a `need_simple` backend**,
  because `std.testing.fuzz` returns immediately on those and every target would
  report green while testing nothing.

The one throughput datum available is an anchor, not a budget: 1,000,007
iterations of a *trivial* target in 15.1 s on aarch64-macos. Coordinate
generation is orders of magnitude heavier per iteration, so per-target limits
cannot be extrapolated from it and must be measured when the harness exists.

## Performance gate — baseline infrastructure, threshold deliberately unset

`zig build performance-baseline -Denable-oracle=true` now records absolute
oracle median and p95 generation time for 20 versioned members in each atom-count
bucket. The harness is fixed to ReleaseFast, runs both warmup and measured calls
in one process, times the public C generation call and result deallocation, and
uses a monotonic clock. CI captures the report on both supported native targets.
Its 100-member population is listed in `tests/oracle_benchmark.zig`: all seven
drug-like members plus the lowest-index adversarial members needed to fill every
bucket. Changing that population requires a decision Bead.

There is still no native implementation to time. Any native/oracle ratio stated
today — including the epic's illustrative 1.5× — would be invented, so none is
set here. `zig build performance-check` therefore remains an explicit failure,
not a false-green gate.

What is fixed now is the **method and the trigger**:

- **Method**: native and oracle timed on the same corpus, in the same build,
  same target, same optimize mode, in one process — the same identity rule
  `RunIdentity.requireSameBuild()` already enforces for coordinates. Report
  median and p95 per size bucket, so a regression confined to `huge` members is
  visible rather than averaged away.
- **Trigger**: the threshold is fixed **when the first representative
  native-vs-oracle baseline exists** — the same milestone that recalibrates T —
  and not later. It is recorded in-repo as a ratio per size bucket.
- **Until then** `cgz-7v2.4.7` records absolute timings each run so the baseline
  accumulates. A missing threshold is not a passing gate: the performance gate
  is `addFail` with its owning bead until the number is set.

Rationale for holding: CoordGen runs interactively inside RDKit depiction, so a
10× regression is a real product failure. But a fabricated threshold is worse
than an unset one — it would either be so loose it never fires, or so tight it
gets relaxed on first contact, which is the `cgz-r21` failure mode.
