# Compatibility and defect policy

Status: frozen by `cgz-r13` against upstream commit
`d20e735d96480385b2e257522288004038a08cc9`. This document makes the epic's
compatibility-and-defect policy reviewable: it defines the four categories a
decision may select, and it defines the fourth one — **parity ceiling** —
narrowly enough that it cannot be used to excuse work.

Observed behavior at the pin is the initial target. When behavior appears
defective, a decision bead selects exactly one category and records the
evidence.

| Category | Meaning | Obligation |
|---|---|---|
| **Emulate** | Upstream's behavior is observable and reproducible; the port reproduces it, defect and all | Oracle regression pinning the behavior |
| **Correct** | Upstream's behavior is a safety, memory, or concurrency defect with no observable single-threaded difference | Evidence that the observable behavior is unchanged, plus a regression for the corrected property |
| **Option** | Consumers plausibly depend on both behaviors | Both paths tested; the upstream-compatible path lives in a named compatibility facade, never in the default |
| **Parity ceiling** | The oracle is not reproducible against itself, so no exact or tolerant claim can be made against it | Enumerated per (member, observable) in the published manifest with the measurement that placed it there |

Never silently change behavior during a mechanical port. Prefer safe corrected
idiomatic behavior, isolating compatibility behavior.

## Parity ceiling

A parity ceiling exists where the **oracle disagrees with itself**. It is not a
category for behavior that is merely hard to reproduce, and it is not a budget.

### The ceiling is the allocator-order axis, not the architecture axis

`cgz-r05` perturbs two variables independently. Only one of them produces a
parity ceiling, and conflating them is the failure mode this section exists to
prevent.

| Axis | What it perturbs | Coordinate divergence (adversarial, n=2000) | Max deviation | Ceiling? |
|---|---|---:|---:|---|
| Heap address order | iteration order of 94 pointer-keyed `std::set`/`std::map` | 1 (0.05%) | 0.580 bond lengths | **Yes** |
| CPU architecture | floating-point evaluation order and precision | 1414 (70.7%) | 49.414 bond lengths | No |

The architecture axis is **not** a parity ceiling. `COMPARISON_SEMANTICS.md`
already requires the differential runner to execute oracle and native in the
same process, the same build, and the same target, and
`RunIdentity.requireSameBuild()` enforces it. Under that rule both sides of
every comparison are the same architecture, so architectural float divergence
is not a difference the runner can ever observe. Its role in the classification
is to prove that cross-build coordinate goldens are forbidden — not to excuse
the port from comparing coordinates.

Reading the 1414/2000 architecture figure as the ceiling would demote 70.7% of
the adversarial corpus to the invariant tier and leave the port proving almost
nothing about the observable it exists to produce. The measurement does not
force that and does not support it.

### The ceiling is an enumerated set, never a fraction

`conformance/parity_expectations.tsv` and the parity ceiling answer different
questions and must not be substituted for one another.

- `parity_expectations.tsv` bounds **how unstable the oracle may be** under each
  perturbation, as a per-partition fraction with deliberate host headroom. It
  is a property of the oracle.
- The **parity ceiling** is the set of (member, observable) pairs the *port* is
  excused from exact or tolerant comparison on. It is enumerated in
  `conformance/parity_ceiling.tsv`, keyed by corpus member identity.

`cgz-r26` moved that enumeration out of the manifest's `order_unstable` column
and into its own artifact. The decision is unchanged; only its home is. The
column was the measurement that produced the rows and stays published evidence,
but the manifest is a per-(architecture, toolchain, optimize-mode) file, so a
gate reading it would change meaning per host — and this section's whole point
is a claim that holds everywhere. `parity_ceiling.tsv` is keyed by
`{partition}/{index}`, which is portable because the corpus generator is a pure
integer function of exactly those two values, and each row pins the SHA-256 of
the member it names so the identity cannot quietly repoint. The consumer
contract is in [COMPARISON_SEMANTICS.md](COMPARISON_SEMANTICS.md).

A fraction can never define the ceiling. A ceiling of `0.002` on 2000 members
permits four failures but does not say which four, so a port could trade a
genuine defect against the allowance and stay green. The excused set is
therefore enumerated, and a pair not in it is a failure.

At the pin the entire ceiling population is four member × observable groups:

| Member | Order-unstable observables | Coordinate deviation |
|---|---|---:|
| `adversarial/917` | `component_transforms`, `component_transforms_set` | 0.000 |
| `adversarial/1538` | `component_transforms`, `component_transforms_set` | 0.000 |
| `adversarial/1588` | `component_transforms`, `component_transforms_set` | 0.000 |
| `adversarial/1695` | `coordinates`, `component_transforms`, `component_transforms_set` | 0.580 |

Nothing else in either partition diverged under heap address order. All
structural observables — ring sets, fragment trees, canonical and input
mappings, morgan ranks, template mappings, effective bond orders, component
membership, `ok`, `clean_pose` — were stable on all 2007 members on both axes,
and remain exact requirements for ceiling members too. The ceiling is scoped to
the named observable, never to the member as a whole.

### What the port must prove about ceiling inputs

Tier 3 (invariant/statistical) from `COMPARISON_SEMANTICS.md` is the **floor**
of the obligation, not the whole of it. Where the manifest publishes a measured
bound on the oracle's own order-instability, that bound is an additional and
tighter obligation:

1. Every observable not named in the member's `order_unstable` cell is compared
   at its normal tier — exact for structural observables, tolerant for stable
   coordinates and derived floats.
2. A named observable is compared tolerantly against the published per-member
   bound where that bound is inside the tolerant band. `adversarial/1695`'s
   0.580 bond lengths is; a port that relayouts that member is a failure, not a
   ceiling hit.
3. Only a pair whose measured order-instability exceeds the tolerant band drops
   to tier 3. **At the pin no pair does, so the invariant/statistical tier is
   currently unpopulated.** Populating it requires a new measurement recorded
   here.

Widening any of this is a weakening under `cgz-r21` and needs a decision bead
naming the measurement that forces it.

### Completion claims

Any completion or parity claim must state this category explicitly and cite the
enumerated population, rather than implying universal exact parity. The correct
form is "exact on all structural observables for all 2007 corpus members; four
member × observable pairs are at the pointer-order parity ceiling and are
enumerated in the parity manifest", not "behavioral parity achieved".

### Superseded figure

`cgz-r13`'s original text, and the epic's summary of it, recorded the
order-axis coordinate divergence as *739 units, roughly 15 bond lengths*. That
came from the standalone C++ review probe in `cgz-r05`. The reproduction
through the Zig corpus generator, the C ABI oracle adapter, and the classifier
measures the same single member at **0.580 bond lengths (29 units)** — a factor
of 25 smaller. The published manifest is authoritative. The corrected figure
strengthens the strict reading: the one input the port structurally cannot
match is off by half a bond length, not by a relayout.

## Register of decided items

| Item | Category | Evidence |
|---|---|---|
| `cgz-r10` — `buildTuplesOfDofs` and upstream 6212c86 | No behavioral delta; nothing to emulate or correct | Upstream issue #137 is a GCC 15 build failure; `runLocalSearch` is uncalled at the pin |
| `cgz-r11` — `initialize()` mutates and takes ownership of caller input | Correct, with the effect surfaced as output | `clear()` at the head of `initialize()` `delete`s the previous call's atoms and bonds |
| `cgz-r12` — template static-init data race | Correct, do not emulate | Non-atomic guard over a mutable static; `normalizeTemplate` is not idempotent |
| `cgz-r13` — pointer-order-unstable inputs | Parity ceiling, scoped as above | 4 member × observable pairs out of 2007, max 0.580 bond lengths |

### cgz-r10 — upstream 6212c86 is not a behavioral change

Commit `6212c86` ("fixes the issue, no idea why") changes one line of
`CoordgenMinimizer::buildTuplesOfDofs` from `for (auto lastOrderTuple : ...)` to
`for (const auto& lastOrderTuple : ...)`. The pin already contains it
(`CoordgenMinimizer.cpp:1089`). Three independent findings settle it:

- **Upstream issue #137 is titled "Doesn't compile with GCC 15."** It reports a
  `memcpy` restriction diagnostic raised inside libstdc++'s
  `uninitialized_copy` while copying a `std::vector`, promoted to an error. PR
  #138 removes the copy and so removes the instantiation. The symptom was
  neither wrong output nor performance; it was a build failure, and the author's
  own words on the PR are that it may be a compiler bug.
- **The function is unreachable at the pin.** Its only caller is
  `CoordgenMinimizer::runLocalSearch`, which has no callers anywhere in the
  pinned tree — not in the library, not in `test/`, not in `example_dir/`. The
  header itself records that `runSearch` is "alternative and preferred to
  runExhaustiveSearch and runLocalSearch". No coordinate-generation path can
  reach the changed line.
- **The change is inert in standard C++.** `lastOrderVector` is not aliased by
  `growingVector` (it is copy-assigned, then `growingVector` is cleared), is not
  mutated in the loop, and the body copies the element into `newTuple`
  regardless. By-value and by-reference iteration produce identical sequences.

`buildTuplesOfDofs` also compares degrees of freedom with pointer **equality**
(`dof == *(lastOrderTuple.rbegin())`), never pointer **ordering**, so its output
order is a pure function of the input vector's order. That is consistent with
the measurement: `dofs`, `dofs_set`, and `dof_penalties` are 0/2000 divergent on
the allocator-order axis. This function is not part of the parity ceiling.

The port does not implement `runLocalSearch`, `buildTuplesOfDofs`,
`runExhaustiveSearch`, or `runExhaustiveSearchLevel` unless a consumer
requirement appears for them, and it carries no analogue of the GCC 15
workaround because it performs no `std::vector` copy.

### cgz-r11 — `initialize()` ownership, and why it is worse than a bond rewrite

`sketcherMinimizer::initialize()` does not merely rewrite bond orders on
caller-owned objects. At the pin it:

1. calls `clear()` as its **first statement**, which `delete`s every atom and
   bond held in `m_referenceAtoms` / `m_referenceBonds` from the previous call
   (`sketcherMinimizer.cpp:338-350`);
2. rewrites `_bond->bondOrder = 0` on caller bonds for non-terminal bonds to
   metals when `getTreatNonterminalBondsToMetalAsZOBs()` is on — **default
   true** (`sketcherMinimizer.h:289`);
3. **erases** zero-order, skipped, and hidden bonds from `minMol->_bonds`, and
   hidden atoms from `minMol->_atoms`;
4. **reorders** both vectors in place via `canonicalOrdering`, which also writes
   `_generalUseN` on caller atoms and `_SSSRVisited` on caller bonds;
5. rewrites the molecule through `forceUpdateStruct` and hands it to
   `splitIntoMolecules`.

Steps 2 and 3 compose: a metal bond demoted to zero order is then deleted from
the caller's molecule outright.

The bead's proposed regression — run one molecule through `initialize()` twice
and diff — **cannot be run**, and that is the finding. The second call's
`clear()` frees the molecule's atoms and bonds, then `m_referenceAtoms =
minMol->_atoms` reads the freed pointers back out and the rest of the function
dereferences them. It is a use-after-free, not a diffable observation.
`initialize()` is therefore not merely non-idempotent: it **takes ownership** of
the caller's atoms and bonds.

**Decision.** The port's public contract is that `generate` never mutates and
never takes ownership of caller input. This is the ownership rule already frozen
in `FOUNDATION_CONTRACTS.md` — Zig input slices and strings are borrowed only
for `generate`, the C implementation borrows every input pointer only for
`coordgen_generate` — and this finding is the concrete defect it protects
against. The metal/ZOB normalization still happens internally and still affects
layout; its result is published as the `effective_bond_order` result array in
caller bond order. Caller-visible reordering and deletion do not happen at all,
because all output arrays are in caller input order by contract.

No compatibility facade reproduces the mutation. A facade could plausibly be
argued for a bond-order rewrite; it cannot be argued for freeing the caller's
objects, and the two are not separable in upstream's implementation.

The replacement regression is: submit the same borrowed input twice, assert the
second result equals the first bit for bit, and assert every caller-owned input
array is unchanged after both calls. That is a stronger claim than the diff the
bead asked for, and unlike it, it is defined behavior.

### cgz-r12 — the template race is real; correct it

Confirmed against the pin, and the mechanism is worse than a double
initialization.

`sketcherMinimizer::loadTemplates()` (`sketcherMinimizer.cpp:3590`) guards
`static CoordgenTemplates m_templates` (`sketcherMinimizer.h:284`, defined at
`sketcherMinimizer.cpp:3750`) with `static bool loaded = false`. Because that
flag is a scalar with constant initialization it is **static-initialized**, so
no `__cxa_guard_acquire` protects it; the read and the write are plain
non-atomic accesses. `findTemplate` calls it from
`CoordgenFragmentBuilder.cpp:110` on every fragment-building path.

Two threads reaching first use together can each observe `loaded == false` and a
non-empty check that passes, and then:

- both execute `m_templates.getTemplates() = coordgen_templates()`, so one
  thread's molecules are leaked and the other thread's in-flight
  `normalizeTemplate` loop walks a vector being reassigned underneath it — a
  use-after-free; or
- both run `normalizeTemplate` over the same molecules. That function is **not
  idempotent**: it divides every atom coordinate by the modal bond length
  (`sketcherMinimizer.cpp:3552-3556`). Running it twice scales templates by
  `1/f²`, producing silently wrong template geometry rather than a crash.

Teardown compounds it: `~CoordgenTemplates` raw-`delete`s every atom, bond, and
molecule at static-destruction time, in an order unspecified relative to other
translation units' statics.

After initialization the data is effectively read-only — `findTemplate` and
`sketcherMinimizer::compare` only read `coordinates` and `_generalUseN` on
template atoms, and `morganScores` writes nothing to them. So making the data
*actually* immutable costs nothing observable single-threaded, which is the
justification the category requires.

**Decision confirmed: correct, do not emulate.** Templates are immutable
compiled data, already normalized at generation time, with no mutable global
initialization and no static destruction order dependency — the requirement
already stated in `FOUNDATION_CONTRACTS.md`. Per-generation state that upstream
stores on template atoms (`_generalUseN`) belongs to the caller's context, not
to the shared data.

The regression is the template-concurrency gate the epic lists: two contexts
generating coordinates on two threads with first use racing, results identical
to the serial case. It is owned by `cgz-7v2.4.3` and cannot be written until the
native port exists.
