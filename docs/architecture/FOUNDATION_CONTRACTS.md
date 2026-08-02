# CoordGen Zig foundation contracts

Status: frozen for the initial behavioral port against upstream commit
`d20e735d96480385b2e257522288004038a08cc9`. Representation changes require a
review bead because topology, layout, optimization, conformance, and interop
all consume these definitions.

## Public execution and ownership

The public API is single-shot: an immutable graph, all stereo, all auxiliary
relations, and all options are submitted together. There is no public sequence
of “add bonds, finalize neighbors, then add stereo”, so the upstream illegal
call order cannot be expressed.

Zig input slices and strings are borrowed only for `generate`. The Zig result
owns its arrays through the caller-supplied allocator and requires exactly one
`Result.deinit`. The C implementation borrows every input pointer only for
`coordgen_generate`; on success the returned `coordgen_result_t` requires
exactly one `coordgen_result_free`. On failure the C result is zeroed. A result
must not be copied and freed twice.

All output atom arrays are in caller input order. All output bond arrays are in
caller input bond order. `input_to_internal` and `internal_to_input` are the
only arrays whose names denote a permutation. Canonical ordering, component
splitting, and fragment construction must never leak their storage order into
consumer output. The adversarial permutation contract is exercised in
`src/model.zig`.

The public bond-length unit is `bond_length = 50.0f` in Zig and
`COORDGEN_BOND_LENGTH = 50.0f` in C.

### Sufficiency evidence

These are claims about what a real consumer needs, so they are settled by a
consumer rather than by argument. `tests/rdkit_consumer.cpp` reproduces the
flow of RDKit's `External/CoordGen/CoordGen.h` — per-atom flags and template
coordinates, bond stereo, one generate call, coordinates read back by the
caller's own atom index and divided by 50 — links it against the pinned C++
facade, and runs it under `zig build upstream-oracle -Denable-oracle=true`.

The molecule is handed over in an order that is deliberately not its
connectivity order, so the input-order contract is falsifiable rather than
assumed: a result returned in canonical or component order would fail the
bonded-distance checks instead of silently transposing coordinates. Measured
against the pinned facade, bonded pairs come back 49.996–50.001 apart, and
renumbering every atom through an unrelated permutation reproduces the distance
matrix with a deviation of exactly 0.

`tests/abi_cpp_consumer.cpp` keeps the same surface as a compile-only check for
builds that must not depend on the oracle; it is the header's portability gate,
not its behavioral one.

## Conserved type table

| Meaning | Zig representation | Units/default | Invalid state | C ABI | Upstream mapping / ownership |
|---|---|---|---|---|---|
| Public atom/bond position | `u32` `InputIndex` | Caller order | `UINT32_MAX`, out of range | `uint32_t` | Replaces caller pointer identity; borrowed input |
| Internal IDs | distinct `enum(u32)` `AtomId`, `BondId`, `RingId`, `FragmentId`, `MoleculeId`, `ResidueId`, `DofId`, `InteractionId` | Dense context-local ordinal | `.invalid = UINT32_MAX` | Probe serializes raw `uint32_t`; not production DTO identity | Replaces typed C++ pointers; context-owned |
| Coordinates/scoring | `f32`, `Vec2`, `Vec3` | CoordGen units; bond = 50 | NaN/Inf rejected at public boundary | fixed-width `float` POD | `sketcherMinimizerPointF` and upstream float evaluation order |
| Atomic number | `AtomicNumber enum(u32)` | 1...118; carbon default | 0 is internal virtual atom and rejected publicly; >118 rejected | `uint32_t` | `sketcherMinimizerAtom::atomicNumber` |
| Formal charge | `i32` | 0 | no sentinel; full signed range conserved | `int32_t` | `sketcherMinimizerAtom::charge` |
| Bond order | `BondOrder enum(u32)` | zero/single/double/triple = 0/1/2/3 | other values rejected | `uint32_t` constants | `bondOrder`; zero order is first-class |
| Atom stereo | `AtomStereo enum(u32)` plus three input references | unspecified; relative clockwise/counter-clockwise or absolute R/S | reference presence must match mode | fixed `uint32_t` plus indices | `m_chiralityInfo`, `isR`, `hasStereochemistrySet` |
| Bond stereo | `BondStereo enum(u32)` plus two substituent references | unspecified; relative cis/trans or absolute Z/E | reference presence must match mode | fixed `uint32_t` plus indices | `sketcherMinimizerBondStereoInfo`, `isZ`, `m_ignoreZE` |
| Bond display | `BondDisplay enum(u32)` | none | unknown values rejected | fixed `uint32_t` | Losslessly normalizes `hasStereochemistryDisplay`, `isWedge`, `isReversed` into none/solid-or-hashed × direction |
| Atom controls | Zig booleans | false | unknown C flag bits rejected | `uint32_t` flag mask | `fixed`, `constrained`, `hidden`, template-coordinate presence, 3D-coordinate presence |
| Residue | `ResidueInput` with atom, chain bytes, residue number, optional closest ligand atom | empty chain/0/no closest atom | bad atom references | pointer-free descriptor plus borrowed string view | `sketcherMinimizerResidue::{chain,resnum,m_closestLigandAtom}` |
| Residue interaction | endpoints, two additional endpoint index slices, crossing multiplier | empty extra endpoint sets; multiplier 1 | bad indices, negative/nonfinite multiplier | flat descriptor and borrowed spans | `sketcherMinimizerResidueInteraction` fields |
| Effective bond order | result array in input bond order | input order unless normalization changes it | n/a | owned `uint32_t` span | Makes the metal/nonterminal zero-order rewrite observable without mutating caller input |
| Clean pose | `bool` | returned for every successful generation | n/a | `uint32_t` 0/1 | `runGenerateCoordinates()` return value |

Public C structs are ordinary structs, never packed structs or C enums with an
implementation-selected width. `tests/abi_layout.c` and Zig layout tests freeze
size, alignment, and offset. All reserved fields must be zero on input.

## Complete upstream option normalization

| Upstream control | Zig field / normalized input | C ABI field | Default | Contract test |
|---|---|---|---|---|
| constructor `precision` | `options.precision`; named `Precision.quick/standard/best` | `precision` | 1.0; named 0.2/1.0/3.0 | API constants/validation |
| `setScoreResidueInteractions` | `score_residue_interactions` | same | true | options default/layout |
| `setTreatNonterminalBondsToMetalAsZOBs` | `treat_nonterminal_bonds_to_metal_as_zero_order` | same | true | options default; effective-order output contract |
| `setEvenAngles` | `even_angles` | same | false | options default/layout |
| `setSkipMinimization` | `skip_minimization` | same | false | options default/layout |
| `setForceOpenMacrocycles` | `force_open_macrocycles` | same | false | options default/layout |
| `setTemplateFileDir` | optional borrowed `template_directory` | string view `template_directory` | null selects generated built-ins | string-view layout/ownership |
| `loadTemplates` | `load_templates` | same | true | options default/layout; implementation must use immutable/context-local data |
| `constrainAllAtoms` | `constrain_all_atoms` | same | false | options default/layout |
| `constrainAtoms(vector<bool>)` | per-atom `constrained` | `COORDGEN_ATOM_CONSTRAINED` | false | lossless atom DTO test |
| `fixAtoms(vector<bool>)` | per-atom `fixed` | `COORDGEN_ATOM_FIXED` | false | lossless atom DTO test |
| `addExtraBond` | `Input.extra_bonds` | `extra_bonds` span | empty | normal validation and caller-index rules |
| `buildFromFragments` | `build_from_fragments` | same | false | options default/layout; conformance/orchestration control |
| `DEBUG_MINIMIZATION_COORDINATES` / review `DEBUG_COORDS` | `debug_coordinates` | same | false | accepted only by builds that provide the debug sink; otherwise `Unsupported` |

The static utilities audited as public (`canonicalOrdering`, `morganScores`,
`RMSD`, `svd`, `alignmentMatrix`, `compare`, `sameRing`, `getBond`) are future
typed functions over normalized inputs or prepared internal graph states, not
raw pointer utilities. `morganScores` must accept a topology-prepared graph (or
return an error), eliminating the upstream `_generalUseN = i` comment-only
precondition.

## DOF and interaction value contracts

`core/dof.zig` covers every subclass present at the pin. Common state is
`{current, optimal, count, tier}` in fixed-width fields. Variable membership is
an `AtomRange` into the owning collection, never a pointer or per-DOF allocation.

| DOF kind | States | Tier | Variant payload / pinned penalty rule |
|---|---:|---:|---|
| flip fragment | root 1, otherwise 2 | 0 | fragment; constrained flip 1000, chain break 10 |
| change parent bond length | 7 | 2 | fragment; 200 × magnitude step |
| rotate fragment | root 1, otherwise 5 | 3 | fragment; 400 × magnitude step |
| scale atoms | 2 | 4 | pivot + affected range; 50 × affected atom count |
| scale fragment | ringless 1, otherwise 5 | 5 | fragment; 500 × magnitude step |
| invert bond | 2 | 1 | pivot, bound + affected range; inverted state 100 |
| flip ring | 2 | 1 | ring, two pivots, affected range, structural multiplier; 200 × multiplier |

DOFs are ordered by `FragmentId`; the three primary fragment DOFs use pinned
upstream order (flip, change-parent-bond-length, rotate), followed by dynamic
DOFs in deterministic builder discovery order over typed IDs. `DofId` is this
stable ordinal. DOF storage belongs to one generation context’s layout/discrete
search lifetime and resets after discrete optimization; probe snapshots copy
into conformance-result-owned storage before reset. The flat
`conformance/probe_types.DofProbe` and matching non-installed
`conformance/include/coordgen_probe.h` record contain kind, IDs, membership
range, current/optimal/count/tier, exact current penalty, and variant
structural penalty.

`core/interaction.zig` covers the five concrete pinned variants: stretch,
bend, clash, positional constraint, and E/Z constraint. Every payload owns only
typed IDs and scalar values. The deterministic interaction list is owned by the
continuous-optimization phase and reset after that phase; nothing exposes a
vtable or raw atom pointer.

## Module ownership and imports

| Layer | Owns | May import |
|---|---|---|
| core | IDs, chemistry/error representations, f32 vectors, bond length, DOF and interaction values | nothing |
| model | working atom/bond/ring/fragment/residue values and order maps | core |
| geometry | geometry interfaces and algorithms | core |
| topology | adjacency, components, rings, canonicalization, CIP/stereo structure | core, model, geometry |
| layout | fragment construction, templates, rings, macrocycles, placement | core, model, geometry, topology |
| optimize | DOF search, interactions, continuous minimization | core, model, geometry; never concrete layout builders |
| generator | orchestration and phase lifetimes | core, model, geometry, topology, layout, optimize |
| api | safe public Zig conversion, validation, owned result | core, generator |
| c_abi | handwritten stable C conversion/exports | core, api |
| conformance | oracle adapters, probes, comparisons | any production layer; never imported by production |

`src/module_layers.zig` is the machine-readable edge allow-list and proves the
graph is one-way. The build graph must create these as named modules and add
only the listed explicit imports. Tests may inspect internals; production code
must never import conformance/oracle modules. Oracle/probe artifacts are never
installed.

## Phase allocation policy

- Caller DTO memory is immutable and borrowed for the public call only.
- `WorkingGraph` owns normalized atoms, bonds, strings, and both permutation
  maps through the caller’s Zig allocator (or the C result owner).
- Phase arenas are explicit in the generator: topology scratch resets after
  topology, layout/DOF storage after discrete search, and interaction scratch
  after continuous minimization. A phase must copy retained values upward.
- Immutable compiled templates contain no mutable global initialization.
- Public validation paths return errors, including OOM. Assertions are only
  for invariants established by prior validation or phase types.
