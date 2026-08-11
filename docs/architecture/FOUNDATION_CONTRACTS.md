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

Borrowed means read-only. `generate` and `coordgen_generate` never write
through a caller pointer and never take ownership of one, so the same input may
be submitted any number of times and every caller-owned array is unchanged
afterwards. This is the contract upstream violates: `sketcherMinimizer::
initialize()` opens with `clear()`, which `delete`s the previous call's atoms
and bonds, then rewrites bond orders, erases zero-order and hidden elements,
and reorders both vectors in place on the caller's molecule. A second
`initialize()` on the same molecule is a use-after-free rather than a
non-idempotent result. `cgz-r11` records this as corrected, not emulated, and
no compatibility facade reproduces it; the layout effect of the metal/ZOB
normalization is published instead through the effective-bond-order output.
See [COMPATIBILITY_POLICY.md](COMPATIBILITY_POLICY.md).

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
| `buildFromFragments` | `build_from_fragments` | same | false; only 0 accepted, nonzero is `Unsupported` | options default/layout; oracle rejects nonzero (cgz-7v2.8) |
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
| topology | adjacency, components, ring *perception*, canonicalization, CIP/stereo structure | core, model, geometry |
| layout | fragment construction, templates, ring *placement*, macrocycles, coordinate assignment | core, model, geometry, topology |
| optimize | DOF search, interactions, continuous minimization | core, model, geometry; never concrete layout builders |
| generator | orchestration and phase lifetimes | core, model, geometry, topology, layout, optimize |
| api | safe public Zig conversion, validation, owned result | core, generator |
| c_abi | handwritten stable C DTOs, and separately the exports that use them | core, api |
| conformance | oracle adapters, probes, comparisons | any production layer; never imported by production |

`src/module_layers.zig` is the machine-readable edge allow-list and proves the
graph is one-way. The build graph must create these as named modules and add
only the listed explicit imports. Tests may inspect internals; production code
must never import conformance/oracle modules. Oracle/probe artifacts are never
installed.

"Rings" appears in two rows because two different things are being owned:
topology decides *which* cycles exist (the perceived ring set, and each ring's
membership and fusion relationships); layout decides *where* they go (template
selection, macrocycle construction, coordinate assignment). A bead implementing
one must not write the other's files.

The table's granularity is one bead per layer, and at that granularity file
ownership is unambiguous — every layer has a distinct root, and
`tools/check-module-imports` already resolves `src/<layer>/` subdirectories, so
a layer can grow files without ambiguity. Splitting a *single* layer across two
parallel beads is not yet specified: `topology`, `layout`, `optimize`, and
`generator` are still `build_support/empty_module.zig` stubs, and until they
have real content there is no sub-file map to divide. Whoever splits a layer
records the file boundary here first.

A layer is a set of modules, not necessarily one. Every module of a layer gets
that layer's row of the table above; a layer's secondary modules additionally
import its canonical module — the one other layers see — under the layer's own
name. `c_abi` is the only layer split this way, and it is split for a linker
reason:

- `c_abi` (`src/c_abi_types.zig`, the canonical module) is the ABI's data
  types and defines no linker name at all.
- `c_abi_exports` (`src/c_abi/exports.zig`) is the sole Zig definition of
  `coordgen_generate` and `coordgen_result_free`, and `src/coordgen.zig` — the
  installed library's root — is the only file permitted to import it.

`conformance/oracle_adapter.cpp` defines those same two symbol names on
purpose, so that a C or C++ consumer of `include/coordgen_abi.h` can link the
pinned upstream oracle without being rewritten. Any binary wanting the
oracle's implementation *and* the ABI's Zig type declarations — the corpus
runner is exactly that — would otherwise link two implementations of the same
ABI. That is a hard `ld.lld: duplicate symbol` error on Linux and silently
tolerated by macOS's linker, so the constraint is enforced statically instead:
`tools/check-module-imports` asserts `c_abi_exports` has exactly one importer,
and the reference is a `comptime` one, because a plain `pub const` `@import`
does not force export discovery in a non-test build.

## Phase allocation policy

- Caller DTO memory is immutable and borrowed for the public call only.
- `WorkingGraph` owns normalized atoms, bonds, strings, and both permutation
  maps through the caller’s Zig allocator (or the C result owner).
- Phase arenas are explicit in the generator: topology scratch resets after
  topology, layout/DOF storage after discrete search, and interaction scratch
  after continuous minimization. A phase must copy retained values upward.
- Immutable compiled templates contain no mutable global initialization, no
  lazy load guard, and no static destruction order dependency. Upstream's
  `loadTemplates()` guards a mutable static with a non-atomic `static bool`
  and then mutates the templates in place through a `normalizeTemplate` that
  is not idempotent, so a first-use race either frees templates under another
  thread or scales their coordinates twice. Normalization happens once, at
  template generation time, and per-generation scratch that upstream stores on
  template atoms belongs to the context. `cgz-r12` records this as corrected.
- Public validation paths return errors, including OOM. Assertions are only
  for invariants established by prior validation or phase types.

## What each gate proves

cgz-r28 was a gate certifying an empty archive: `check-install-isolation`
asserted that `libcoordgen.a` existed, that the probe header did not, and that
no upstream C++ symbols had leaked — all true of a 2400-byte file exporting
nothing. It never asked whether the library contained the ABI it claims to
install. That is the third instance of one pattern (cgz-r19 was flags weakened
to pass, cgz-r20 was steps reporting success unconditionally), so the gates are
enumerated here with the distinction each one turns on.

**Behavioral — recomputes or executes the thing it certifies:**

| Gate | Proves |
|---|---|
| `tools/verify-upstream` | The pinned archive's SHA-256, byte count, Git tree hash, LICENSE hash, and Zig package hash, all recomputed with the pinned toolchain. The strongest gate in the tree. |
| `tools/check-install-isolation` | Builds into a scratch prefix created empty and destroyed after, then compiles and runs `tests/install_consumer.c` against `$prefix` alone. An archive without the ABI fails to link. Wired into `package-check`. |
| `abi-check` | Links `tests/abi_layout.c` against the real library, so a name or signature drift from `include/coordgen_abi.h` is a link error. |
| `assertNoFalseGreenSteps` (build.zig) | Asks the constructed graph whether each public step has dependencies. Runs at configure time on every invocation. |
| `tools/check-gate-strength` | Every one of build.zig's C/C++ compile sites uses an approved strict flag set, with a narrow enumerated exemption for upstream-owned sources. Universality, not presence. |
| `tools/check-module-imports` | No relative import crosses a layer; no bare import lacks an approved edge; `c_abi_exports` has exactly one importer. |

**Textual — reads the build description, and only proves what it says:**

`tools/check-build-policy` asserts that build.zig and build.zig.zon *declare*
a lazy oracle dependency, an `enable-oracle` opt-in, a public-header install,
and no probe-header install. Those are properties of the build description, so
text is the right level for them — but a declaration is not an outcome, and
its summary line says so. The behavioral counterparts are named above.

Both checkers with a `--self-test` mode (`check-module-imports`,
`check-gate-strength`) plant a violation of each class they claim to catch,
plus a compliant control, so a gate cannot decay into one that passes
everything. New gates follow the same rule: state what is asserted, and if the
name implies content, assert content.
