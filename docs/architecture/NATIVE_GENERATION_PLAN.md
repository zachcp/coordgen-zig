# Native coordinate-generation translation plan

Status: implementation decomposition for upstream commit
`d20e735d96480385b2e257522288004038a08cc9`. This plan is owned by
`cgz-7v2.11`. The representation, ownership, comparison, and compatibility
contracts in the existing architecture documents remain authoritative.

## Why this plan exists

The Phase 1 geometry and graph foundations landed as `cgz-7v2.9` and
`cgz-7v2.10`, but the Phase 2–5 implementation beads promised by `cgz-7v2`
were never created. When this document was written the `layout`, `optimize`,
and `generator` build modules were empty-module stubs. This document closes
that planning gap. It maps the pinned C++ implementation to one native owner
per component, defines the integration order, and names the conformance
evidence required at each seam.

## Current state (2026-08-16)

Every module this plan decomposes is implemented and merged. `cgz-7v2.12`
through `cgz-7v2.20` all landed in PRs #19 and #20 and are closed; `src/` is
~19,400 lines and CI is green on `main`.

Integration landed with `cgz-7v2.21`: `api.generate` runs the full pipeline
and `coordgen_generate` publishes owned spans for all six observables, so a
supported input now returns coordinates through both the safe and the C entry
point. Coverage is partial **by domain, never by size** — the rejected set is
enumerated on `api.generate` and an unowned domain returns `Unsupported`
rather than a successful empty result.

What has still not happened is comparison. **No native output has yet been
compared against the oracle** — the tolerance `T` in
[`SUCCESS_CRITERIA.md`](SUCCESS_CRITERIA.md) is still a prior taken from the
oracle's own float sensitivity, and the performance gate is still deliberately
red with no threshold.

The beads for this decomposition were created retroactively on 2026-08-16,
eleven days after the work merged. See `cgz-r29` for why that happened and
what should prevent the next one.

In the tables below, `P/` means the pinned package root recorded by
`upstream/coordgenlibs.lock` and materialized by Zig under its content hash.
Line ranges name the audited pin; function names are the durable mapping when
generated or surrounding lines move.

## Upstream execution order

### Input preparation

`sketcherMinimizer::initialize` (`P/sketcherMinimizer.cpp:164-280`) performs
these observable steps:

1. Clear prior context and retain caller order.
2. Rewrite eligible nonterminal metal bonds to effective zero order.
3. Convert zero-order, extra, and residue-interaction bonds into proximity
   relations.
4. Exclude skipped, zero-order, and hidden graph elements from structural
   adjacency.
5. Compute Morgan/canonical order.
6. Populate working atoms, bonds, residues, and residue interactions.
7. Assign neighbors and perceive rings.
8. Split connected components.
9. Attach proximity relations to components and flag cross-layout atoms.
10. Publish the prepared graph to the minimizer.

The native implementation preserves the effects but not upstream's destructive
ownership. It copies borrowed input into `model.WorkingGraph`, records effective
bond order as output, and keeps both caller/internal order maps as required by
`cgz-r11` and `FOUNDATION_CONTRACTS.md`.

### Generation

`sketcherMinimizer::runGenerateCoordinates`
(`P/sketcherMinimizer.cpp:295-312`) is the reference phase order:

```text
structurePassSanityCheck
  │
  ├─ findFragments
  │    ├─ CoordgenFragmenter::splitIntoFragments(component)
  │    └─ initializeFragments
  │         ├─ assignNumberOfChildrenAtomsFromHere
  │         ├─ CoordgenFragmentBuilder::initializeCoordinates(fragment)
  │         └─ assignLongestChainFromHere
  │
  ├─ buildFromFragments(true)
  │    └─ CoordgenMinimizer::buildMoleculeFromFragments
  │         ├─ alignWithParentDirection
  │         └─ fragment.setCoordinates
  │
  ├─ CoordgenMinimizer::avoidClashes
  │    ├─ build and score clash/constraint interactions
  │    ├─ flipFragments / tiered discrete DOF search
  │    ├─ rebuild fragment coordinates
  │    ├─ avoid terminal clashes
  │    └─ continuous minimization when required
  │
  ├─ bestRotation
  ├─ maybeFlip
  ├─ arrangeMultipleMolecules
  ├─ writeStereoChemistry
  └─ return aggregate clean_pose
```

Native `generator` owns this sequence and explicit phase lifetimes. A phase
copies retained state to its caller before resetting its scratch arena.

## Component coverage and ownership

Each production method belongs to exactly one row. Public declarations in C++
headers describe the same implementation and do not create a second owner.

| Component | Pinned upstream implementation | Native owner and normalized result | Conformance evidence |
|---|---|---|---|
| Input normalization and working graph | `sketcherMinimizer.cpp:65-280,314-506,3609-3651` (`canonicalOrdering`, `initialize`, `flagCrossAtoms`, `splitIntoMolecules`, `sameRing`, `morganScores`); `sketcherMinimizerMolecule.{h,cpp}` all; atom construction/basic valence/ring/coordinate/submolecule helpers in `sketcherMinimizerAtom.cpp:204-373,790-806,1220-1235,1506-1680` | `model` and `topology`: immutable-input conversion, effective orders, adjacency, components, ring perception, canonical maps | Exact input/internal maps, Morgan ranks, component spans, ring memberships; hidden/skipped/zero-order/metal and order-permutation properties |
| Topology and stereochemistry | CIP helper and atom stereo families in `sketcherMinimizerAtom.cpp:21-203,351-789,807-1219,1236-1505`; `sketcherMinimizerBond.{h,cpp}` all; `sketcherMinimizerRing.{h,cpp}` all; `sketcherMinimizer.cpp:507-611` (`writeStereoChemistry`, `assignPseudoZ`, peptide stereo flips) | `topology`: ring chemistry/fusion, CIP traversal and priority, absolute atom/E/Z state, final display descriptors. Ring *placement* is excluded. | Exact ring structure through probes; atom stereo, bond stereo, and display through the stable result ABI; mirrored-chiral negative fixture |
| Fragmentation and basic placement | `CoordgenFragmenter.{h,cpp}` all; `CoordgenFragmentBuilder.cpp:29-102,166-793,813-1090` except template and macrocycle dispatch identified below; fragment container/frame methods `sketcherMinimizerFragment.cpp:366-525`; `sketcherMinimizer.cpp:2715-3055`; shared coordinate helpers in `sketcherMinimizerMaths.h` | `layout`: fragment membership/priority/parent tree, acyclic placement, regular and fused-ring placement, constrained parent alignment, local coordinate frames | Exact fragment atoms/rings/parent/component/flags; tolerant coordinates; upstream basic, bicyclic, fused-ring, and fragment-membership cases |
| Templates | `CoordgenTemplates.{h,cpp}` all (generated production data); `CoordgenFragmentBuilder::findTemplate` at `cpp:103-165`; constraint setters plus RMSD/alignment/SVD/compare/identity/load families at `sketcherMinimizer.cpp:3056-3608` | `layout/templates`: immutable generated template data, graph identity/mapping, alignment, and coordinate application. No mutable first-use state. | Exact template index and atom mapping; all 82 normalized templates; byte-identical regeneration; tolerant resulting coordinates; concurrent first use equals serial |
| Macrocycles | `CoordgenMacrocycleBuilder.{h,cpp}` all; `CoordgenFragmentBuilder::maybeAddMacrocycleDOF` at `cpp:794-812`; macrocycle branches of `buildRing` and `buildRings` | `layout/macrocycle`: Polyomino/Hex state, path/ring/EZ constraints, shape enumeration/scoring, chosen open bond/path, coordinates and ring DOF | Exact ring/fragment and flip-ring DOF records; Polyomino public tests; `macrocycle.mae`; tolerant coordinates and invariant quality where oracle choice is unstable |
| Discrete optimization | All seven DOF implementations in `sketcherMinimizerFragment.cpp:29-365`; `CoordgenMinimizer.cpp:719-1433,1638-1695` (scores, exhaustive/tuple/local/tiered search, solutions, flips, internal/terminal clash handling) | `optimize/discrete`: value-owned DOF application, deterministic tuple/search order, clash/cross-bond/ring/proximity scores, selected state | Exact DOF kind, IDs, states, tier, affected atoms and pivots; tolerant penalties; exact `clean_pose`; all-index OOM and deterministic search properties |
| Continuous optimization | `CoordgenMinimizer.cpp:43-718,1491-1637`; all six `sketcherMinimizer*Interaction.h` implementations; `sketcherMinimizerMarchingSquares.{h,cpp}` all | `optimize/continuous`: interaction construction, scoring/forces, precision/termination, NaN and 3D fallback, residue/molecule minimization, scalar contours. Marching Squares has this single owner even though it is off the minimal path. | Tolerant final coordinates/energies where observable; stretch/bend fixtures; malformed/nonfinite and OOM properties; corpus quality invariants |
| Orchestration and finalization | Remaining controller code in `sketcherMinimizer.{h,cpp}`: constructor/config/sanity/run `cpp:48-64,282-337,3691-3750`; global orientation `612-1020`; residues, chains/SSEs, crowns, proximity grids, molecule placement and crossing flips `1074-2713`; `sketcherMinimizerResidue.{h,cpp}` and `sketcherMinimizerResidueInteraction.{h,cpp}` all | `generator`: phase sequencing/lifetimes, sanity/fallback, global orientation, residues/proximity, multicomponent placement, final stereo invocation and clean-pose aggregation | Exact component memberships/transform status/clean pose; tolerant component transforms; residue and disconnected corpus partitions; repeat/thread/order determinism |
| Safe and C APIs | Upstream public declarations in model/controller headers are interface reference only. There is no upstream C ABI. | `api`, `c_abi`, and installed library root: borrowed DTO conversion, allocator-owned result, caller-order serialization, errors, `coordgen_generate` and `coordgen_result_free` | Safe API tests; strict ABI/layout/errors; successful allocation/free and injected OOM; installed Zig/C/C++/RDKit-shaped generation consumers |

### Pinned top-level source ledger

This ledger makes coverage mechanically auditable by exact basename. A row owns
both declarations and definitions unless the notes explicitly partition a
controller file by the ranges in the component table above.

| Exact pinned source | Classification / single translation owner |
|---|---|
| `CoordgenFragmentBuilder.cpp`, `CoordgenFragmentBuilder.h` | Fragment/basic layout, except the exact template and macrocycle method ranges assigned above |
| `CoordgenFragmenter.cpp`, `CoordgenFragmenter.h` | Fragment/basic layout |
| `CoordgenMacrocycleBuilder.cpp`, `CoordgenMacrocycleBuilder.h` | Macrocycles |
| `CoordgenMinimizer.cpp`, `CoordgenMinimizer.h` | Continuous or discrete optimization by the non-overlapping function ranges above; fragment assembly entry points belong to basic layout |
| `CoordgenTemplates.cpp`, `CoordgenTemplates.h` | Templates; generated production data plus its declarations |
| `sketcherMinimizer.cpp`, `sketcherMinimizer.h` | Partitioned controller: normalization, stereo, layout alignment, templates, or orchestration by the exact ranges above; no unassigned production range |
| `sketcherMinimizerAtom.cpp`, `sketcherMinimizerAtom.h` | Normalization/basic graph helpers or topology/stereo by the exact ranges above |
| `sketcherMinimizerBond.cpp`, `sketcherMinimizerBond.h` | Topology/stereo |
| `sketcherMinimizerRing.cpp`, `sketcherMinimizerRing.h` | Topology/stereo structure; layout consumes but does not reimplement it |
| `sketcherMinimizerFragment.cpp`, `sketcherMinimizerFragment.h` | Discrete optimization owns DOF methods; fragment/basic layout owns container and coordinate-frame methods |
| `sketcherMinimizerInteraction.h`, `sketcherMinimizerBendInteraction.h`, `sketcherMinimizerClashInteraction.h`, `sketcherMinimizerConstraintInteraction.h`, `sketcherMinimizerEZConstrainInteraction.h`, `sketcherMinimizerStretchInteraction.h` | Continuous optimization |
| `sketcherMinimizerMarchingSquares.cpp`, `sketcherMinimizerMarchingSquares.h` | Continuous optimization |
| `sketcherMinimizerMaths.h` | Fragment/basic layout; other components consume the frozen geometry seam |
| `sketcherMinimizerMolecule.cpp`, `sketcherMinimizerMolecule.h` | Input normalization/working graph |
| `sketcherMinimizerResidue.cpp`, `sketcherMinimizerResidue.h` | Orchestration/finalization |
| `sketcherMinimizerResidueInteraction.cpp`, `sketcherMinimizerResidueInteraction.h` | Orchestration/finalization |
| `sketcherMaeReading.h` | Parsing-only; explicitly outside native generation |

### Classified outside the native generator

- `sketcherMaeReading.h` is parsing-only reference code; native generation
  accepts normalized DTOs and does not own MAE parsing.
- `test/*` and `example_dir/example.cpp` are tests/examples mapped by the
  requirements-coverage manifest, not production translation units.
- `CoordgenTemplates.cpp` is generated production data and belongs to the
  template component rather than the handwritten algorithm inventory.
- `sketcherMinimizer::writeMinimizationData` is debug-only. A native debug sink
  remains opt-in and must not alter production results.
- `CoordgenMinimizer::runLocalSearch` is dormant in the pinned main path but is
  production-capable and remains owned by discrete optimization.

## Native module boundaries

The canonical roots remain the modules named by `src/module_layers.zig`.
Submodules below each root provide exclusive implementation ownership:

```text
src/model.zig                         long-lived working values and order maps
src/topology.zig                      topology public seam
src/topology/prepare.zig              normalization and effective orders
src/topology/canonical.zig            Morgan/canonical ordering
src/topology/rings.zig                ring perception and fusion structure
src/topology/stereo.zig               CIP, chirality, E/Z, display
src/layout.zig                        layout public seam
src/layout/fragments.zig              fragment graph and priority/tree
src/layout/basic.zig                  acyclic and ordinary ring placement
src/layout/templates.zig              immutable templates and matching
src/layout/macrocycle.zig             Polyomino search and macrocycle placement
src/optimize.zig                      optimization public seam
src/optimize/discrete.zig             DOF application/search and clash choices
src/optimize/continuous.zig           interactions, forces, minimization
src/optimize/marching_squares.zig     scalar contour helper
src/generator.zig                     orchestration public seam
src/generator/context.zig             phase storage and reset boundaries
src/generator/residues.zig            residue/proximity placement
src/generator/components.zig          global orientation/component placement
src/api.zig                           safe public conversion and owned result
src/c_abi/exports.zig                 stable ABI conversion/ownership only
```

Only the bead named as owner may create or edit a listed submodule. Canonical
root and build-graph wiring changes belong to the integration bead, after the
submodule implementations settle; this prevents parallel feature beads from
editing the same roots. Existing foundation files may be extended only by the
bead that explicitly owns that extension.

## Implementation sequence

The critical path deliberately reaches a successful, useful native call before
the full upstream feature set:

1. **Prepare graph.** Convert a visible connected acyclic input, derive effective
   bond orders, adjacency, components, and canonical maps.
2. **Basic layout.** Fragment and place an acyclic stereo-free molecule at bond
   length 50 with deterministic local coordinates.
3. **Minimal internal integration.** Assemble one component and map coordinates
   back to caller order through an internal generator seam and executable test.
   This proves real native coordinate generation without prematurely editing the
   canonical roots or promising unsupported public domains.
4. **Topology/stereo and ordinary rings.** Add complete ring and stereo behavior
   without changing the public lifetime contract.
5. **Templates, then macrocycles.** Add immutable matching/placement and the
   specialized cyclic search.
6. **Discrete and continuous optimization.** Implement all DOF choices, clashes,
   interactions, and minimization.
7. **Residues, proximity, multicomponent placement, and final orchestration.**
8. **Close native-dependent gates.** Run differential, OOM/concurrency, installed
   consumers, platform, fuzz, and performance evidence against real generation.

The minimal success fixture is a connected acyclic stereo-free molecule with
no template, residue, proximity, macrocycle, constraint, clash, or minimization
requirement. The internal result must contain deterministic caller-order
coordinates at bond length 50. Public `generate` remained unsupported until the
full integration bead could distinguish every supported domain and provide the
complete ownership contract; unsupported domains must never return a silently
empty successful result. `cgz-7v2.21` discharged that condition: both entry
points reject by domain with `Unsupported`, and the C result publishes all six
spans under one owner released by `coordgen_result_free`.

## Implementation beads and dependencies

The Beads database is authoritative for current IDs and status. The planned
dependency graph is:

```text
cgz-7v2.11  translation plan
    │
    ▼
cgz-7v2.12 prepare graph
    │
    ▼
cgz-7v2.13 topology/stereo
    │
    ▼
cgz-7v2.14 fragment/basic layout
    │
    ▼
cgz-7v2.15 minimal internal generation
    │
    ├───────────────┬──────────────────┐
    ▼               ▼                  ▼
cgz-7v2.16      cgz-7v2.17         cgz-7v2.18
templates       macrocycles        continuous optimize
    │               │                  │
    └───────────────┴──────────┬───────┘
                               ▼
                        cgz-7v2.19
                        discrete optimize
                               │
                               ▼
                        cgz-7v2.20
                  residues/components/finalization
                               │
                               ▼
                        cgz-7v2.21
                    full-generation integration
```

Dependencies encode ordering where a consumer requires settled semantics.
Independent submodules may proceed in parallel because their file scopes do not
overlap. The final integration bead alone edits canonical roots, build wiring,
and public exports after all feature modules are available.

## Required evidence per implementation bead

Every child bead must provide all of the following before closure:

1. Mapped upstream files, function families, and any excluded branch.
2. Normalized inputs, retained outputs, allocator ownership, and scratch reset.
3. Exact intermediate probe comparison wherever the observable is stable.
4. Tolerant coordinate comparison in bond-length units and no reflection by
   default.
5. Invariant/statistical evidence only where the stability catalog justifies a
   parity ceiling.
6. Direct unit/property/malformed/OOM tests for the owned seam in Debug and
   ReleaseSafe from fresh local caches.
7. At least one falsification demonstrating that the gate fails when its owned
   result is perturbed or omitted.
8. No weakening of flags, thresholds, corpus, assertions, or reserved gates.

## Validation unblocking

The full-generation integration bead is the concrete prerequisite that was
previously missing from:

- `cgz-7v2.4.2` differential runner,
- `cgz-7v2.4.3` generation-level property/OOM/concurrency coverage,
- `cgz-7v2.4.4` bounded full-generation fuzzing,
- `cgz-7v2.4.6` installed native generation and ownership,
- `cgz-7v2.4.7` native/oracle performance and quality ratios, and
- `cgz-r16` successful C-entry allocation/free and injected OOM probes.

Those dependencies must point to the concrete integration bead. Their existing
non-generation slices remain valid and must not be reopened.

As of 2026-08-16 that integration bead is `cgz-7v2.21`, it exists, and all six
edges above are recorded in the database. Each of the five validation beads
also carries a note distinguishing the slice that already merged from the slice
that waits on generation, so none of them is reopened wholesale.

## Remaining sequence

1. **`cgz-7v2.21` — wire public generation.** Nothing else can start. The
   ABI/consumer tests that currently assert `unsupported` are part of this
   diff, and unsupported domains must return a specific error rather than a
   silently-empty success.
2. **`cgz-r26` — give the enumerated parity ceiling a portable home.** Wanted
   before the differential runner, not after: the runner has to know per member
   whether an exact claim is owed, and neither existing artifact can tell it
   portably. This is decidable today and does not wait on step 1.
3. **`cgz-7v2.4.2` — differential runner, then the first native baseline.**
   That single baseline discharges two deferred numbers: it recalibrates `T`
   (which today bounds the *oracle's* float sensitivity, not the port's) and it
   sets the per-bucket performance ratio for `cgz-7v2.4.7`. Lowering `T` is
   free; raising it needs a decision bead.
4. **Fan out the remaining validation slices** — `cgz-7v2.4.3` generation-level
   property/OOM/concurrency, `cgz-7v2.4.4` fuzzing (which also owns failure 2
   of `cgz-r27`, the 0-byte crash seed), `cgz-7v2.4.6` installed native
   generation, `cgz-7v2.4.7` ratios, `cgz-r16` C-entry OOM probes. These are
   independent of each other once step 3 exists.
5. **`cgz-7v2.4.8` — final validation audit and scoped parity report.**

`cgz-r25` is settled: the `build_from_fragments` slot is now
`coordgen_options_t.reserved`, required to be zero, and gone from
`api.Options` entirely, so step 1 freezes public generation around a reserved
field rather than around a name that promises a toggle. See
[FOUNDATION_CONTRACTS.md](FOUNDATION_CONTRACTS.md). One open item remains cheap
now and expensive later: `cgz-7v2.3` (owner-controlled mirrors for the pinned
archives, which every gate above assumes will keep resolving).
