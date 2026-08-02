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
