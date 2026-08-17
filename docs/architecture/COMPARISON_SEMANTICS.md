# Differential comparison semantics

`cgz-r06` freezes the policy for native-vs-oracle differential tests. The
executable contract is
[`src/conformance/comparison.zig`](../../src/conformance/comparison.zig); this
document states the policy later runners must preserve.

## Assignment and tiers

The stability catalog assigns a tier per **input × observable**, not per input
or phase. A stable ring record therefore remains exact even when coordinates
for the same input are unstable.

| Tier | Applies to | Required comparison |
|---|---|---|
| Exact | Stable discrete topology, mappings, fragment/DOF state, flags, stereo, displays, templates, components, `ok`, and `clean_pose` | Byte/value equality |
| Tolerant | Stable coordinates and derived floats such as penalties and component transforms | Translation/rotation-normalized geometry or numeric deviation, in bond-length units (`bond_length = 50`) |
| Invariant/statistical | Any input-observable pair the oracle cannot reproduce under catalog perturbations | Corpus-aggregated quality floor: clash score, bond-length RMS, bond-angle deviation, crossings, atoms inside rings, and `ok` |

For the invariant tier, native must be no worse than the same-build oracle plus
a named per-metric margin. Margins are test data, never an implicit global
tolerance.

This document fixes *how* to compare. The numbers — the tolerant-tier tolerance,
the invariant-tier margins, and the pass counts per tier — are set by `cgz-r14`
in [`SUCCESS_CRITERIA.md`](SUCCESS_CRITERIA.md). Note that the tier assignment
inputs (`conformance/parity_expectations.tsv`) bound oracle against *perturbed
oracle*; they are not native-vs-oracle tolerances.

`SUCCESS_CRITERIA.md`'s tier-3 section still reads "on this corpus, coordinates
on the 1414 arch-divergent adversarial members". That sentence predates
`cgz-r13`, which rejected exactly that reading: demoting 70.7% of the
adversarial corpus to the invariant tier would leave the port proving almost
nothing about the observable it exists to produce, and the same-build rule means
the runner never sees the architecture axis at all. `cgz-r13` governs. The
invariant/statistical tier is **unpopulated** at the pin. Correcting `cgz-r14`'s
text is a tightening and needs its own bead.

## The differential runner's per-pair rule

`cgz-r26` settles which artifact the runner reads. There are three, and only
one of them is authoritative for the parity ceiling.

| Artifact | Keyed by | Portable | Authoritative for |
|---|---|---|---|
| `conformance/parity_expectations.tsv` | (partition, observable) | yes | how often the *oracle* may disagree with itself, as a fraction |
| `conformance/parity_manifest.tsv` | corpus member, per build | **no** | published evidence for the build that produced it |
| `conformance/parity_ceiling.tsv` | (member, observable) | yes | **the parity ceiling** |

The expectations file cannot carry the ceiling: its schema is per-partition
fractions and a fraction permits *N* failures without naming which *N*, which is
the budget framing `cgz-r13` rejected. The manifest names members but is a
per-(architecture, toolchain, optimize-mode) artifact, so a gate reading it
changes meaning per host. The ceiling file exists because the enumeration needs
both properties at once: named members *and* the same meaning everywhere.

A member is identified as `{partition}/{index}` — the key
`tests/oracle_corpus_run.zig` and the manifest already use. That identity is
portable because `src/conformance/corpus.zig` generates every member from
nothing but those two values: the xorshift64\* PRNG is pinned by value rather
than taken from a std default, `generateAdversarial` seeds `Rng.init(index)` so
no member depends on any other, and every draw is integer arithmetic and modulo
reduction over a fixed element table. There is no float, no clock, no file, and
no allocator-ordered iteration in generation, so a toolchain bump, an
optimize-mode change, and a different architecture all produce the same
molecule. The one thing the identity does not survive is an edit to the
generator, which would leave every row parsing while silently naming a
different molecule — so each row also carries the SHA-256 of the member's
canonical dump, and `zig build parity-ceiling-check` regenerates the member and
compares. That gate needs no oracle, which is the point.

`comparison.differentialComparison` is the executable form of the rule:

1. Order-stable pair → its ordinary tier. Exact for structural observables,
   tolerant at `T` for coordinates and derived floats.
2. Order-unstable pair with an enumerated row → tolerant, at the row's
   published bound if it has one and at `T` if it does not. A row without a
   bound excuses the pair from *exact* comparison and widens nothing.
3. Order-unstable pair with no row → `error.UnenumeratedOrderInstability`.
   **Fail closed.** Falling through to the invariant tier would let any newly
   unstable pair demote itself to the weakest comparison in the policy, which
   is exactly what enumerating the ceiling prevents.

Stability here means the **heap-address-order axis alone**. The architecture
axis is not an input: the same-build rule below means both sides of every
comparison share a target, so architectural divergence is not a difference the
runner can observe. `tierFor` predates this rule and reads
`Stability.exact()`, which is both axes; runners must use
`differentialComparison`.

The file cannot express a weakening. Its loader rejects an observable whose
normal tier is exact — structural observables stay exact on ceiling members
too — rejects a bound on any observable but `coordinates`, since that is the
only one the classifier measures a per-member magnitude for, and rejects a
bound at or above 1.0 bond lengths, the boundary above which a "tolerant" pass
no longer means the layout matched. Admitting a member is a review event with a
recorded measurement, the same standing as regenerating the manifest.

`zig build corpus-check` enforces the same rule against the measurement:
any (member, observable) pair the classifier finds order-unstable and the
enumeration does not name fails the build by name. A row whose pair is stable
on some host is not a failure — a row is permission, not obligation.

## Geometry and chirality

Coordinate normalization may translate and rotate. It must not reflect by
default: reflection can convert a chemically wrong mirrored layout into a
perfect score. `ReflectionPolicy` permits it only when the fixture records a
non-empty achirality proof. The conformance test intentionally compares a
mirror image and verifies default rejection.

## Baseline identity

Oracle coordinate baselines are artifacts of one `(process, target, toolchain,
optimize-mode)` run; they are not portable golden fixtures. A differential
runner must execute oracle and native in the same process, build, and target,
and `RunIdentity.requireSameBuild()` rejects a mismatch before coordinate
comparison. The stability catalog stays portable because it records category
ceilings instead of raw coordinate goldens.

`runGenerateCoordinates()` / `clean_pose` is a first-class exact output; it is
neither inferred from coordinates nor treated as diagnostic-only.
