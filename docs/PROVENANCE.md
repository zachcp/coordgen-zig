# Upstream semantic provenance

Behavioral compatibility is fixed to Schrödinger CoordGen commit
`d20e735d96480385b2e257522288004038a08cc9`, whose Git tree is
`ec93c5f64e37236cc2e1498732dfdade2eab54dc`. This commit was `master` HEAD
when audited on 2026-08-01, but no build or verification step resolves that
branch name. Changing this commit changes the parity target and requires a
separate review and full re-audit.

`upstream/coordgenlibs.lock` records three independent integrity identities:

- the downloaded archive SHA-256;
- the reconstructed Git tree SHA-1;
- Zig 0.17's normalized package-tree hash used by `build.zig.zon`.

For this URL archive, Zig 0.17 strips the single full-commit top-level
directory from the dependency root. The verifier intentionally hashes that
extracted source root so it reproduces clean-cache URL dependency semantics.
Do not pre-seed Zig's global cache with a differently wrapped directory: the
wrapper changes both the content hash and visible package paths.

`tools/verify-upstream` checks all three, checks the archive root names the
full commit, verifies the preserved license bytes, rejects branch-resolving
URLs, and ensures `build.zig.zon` carries the same URL/hash as a lazy
dependency.

The range `v3.0.2..d20e735` contains ten commits: nine non-merge changes and
the final merge commit. Full SHAs and subjects are enumerated in the lock;
the distinction corrects the ambiguous phrase "3.0.2 plus nine commits."

The upstream BSD-3-Clause license is preserved verbatim at
`LICENSES/coordgenlibs-BSD-3-Clause.txt` and is also required to be present in
the fetched package.

## Template provenance and normalization

The parity target is the normalized template set loaded by upstream, not raw
`templates.mae`. The pinned MAE contains 82 structures, 1704 atoms, and 1963
bonds. `tools/extract-template-reference.py` independently parses the literal
tuples and bonds in committed upstream `CoordgenTemplates.cpp` and produces
`conformance/fixtures/templates_normalized.zig`; it does not use the Zig MAE
reader. A separate raw-mode comparison proves every unnormalized atom tuple
and bond agrees before normalization. The raw and normalized datasets agree
exactly, so no discrepancy preference was required.

Normalization reproduces `sketcherMinimizer.cpp`'s evaluation order: f32
squared bond lengths are promoted for strict double-literal 0.9/1.1 comparisons, the
largest cluster wins with the first cluster winning ties, and all coordinates
are divided by the f32 square root of its representative. Atom indices remain
storage order. The committed representation is immutable Zig data with f32
coordinates encoded as bit patterns, avoiding textual floating-point drift.

`zig build template-regeneration-check -Denable-oracle=true` extracts the C++
reference afresh, generates from `templates.mae` twice, and compares all three
outputs byte-for-byte with the committed fixture. The gate runs natively in
both macOS arm64 and Linux x86_64 CI; fixture equality on each native target
therefore supplies the cross-target determinism proof. Without the explicit
option it fails clearly and the lazy upstream package is not fetched.

## Mirror requirement

The exact archive URL currently points to the full GitHub commit and is
content-checked, but an owner-controlled mirror cannot be named from this
repository: no remote or artifact store is configured. The lock records
`UNCONFIGURED`, verification warns, and `CGZ_REQUIRE_UPSTREAM_MIRROR=1` turns
that condition into an error. Because Zig treats the package hash as the
identity, the URL can later move to the owner-controlled mirror without
changing package contents.
