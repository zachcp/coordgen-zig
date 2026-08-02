# Package and dependency policy

`build.zig.zon` is the only Zig package/dependency manifest. Compiler URLs and
distribution checksums live in the separate host-readable `toolchain.lock`.

The production `coordgen` module starts standard-library-only. Its build
module explicitly forbids libc and libc++, and no dependency module is added
to its import table. Dependency types may never appear in the public API.

The sole initial dependency is the pinned upstream C++ conformance oracle. It
uses an exact normalized content hash and `.lazy = true`. `build.zig` calls the
0.17 `dependencyLazy` API only when `-Denable-oracle=true` is configured.
Clean-cache URL fetching strips the archive's full-commit top-level directory;
the dependency root therefore exposes `LICENSE` and sources directly.
Because the Maker configurator cannot inspect requested top-level step
arguments, the explicit option is the stable boundary:

```sh
zig build                              # never requests the oracle
zig build test                         # never requests the oracle
zig build -Denable-oracle=true conformance
```

The current `upstream-oracle` and `conformance` steps validate package and
license availability; the oracle-build bead will attach compiled artifacts and
tests to these already-stable step names.

## Package contents and offline operation

The `.paths` allowlist includes build files, production sources, public
headers, required examples/assets, README, documentation, provenance locks,
license notices, and build-policy verifiers. Local conformance ABI types,
headers, and compile probes are packaged only because declared build checks use
them; `coordgen_probe.h` is never installed. Upstream C++ conformance sources
remain inside the lazy dependency and are absent from ordinary consumers.

Fetch the complete dependency tree once, then builds can use only local
content-addressed storage:

```sh
zig build --fetch=all
HTTPS_PROXY=http://127.0.0.1:9 zig build test
HTTPS_PROXY=http://127.0.0.1:9 zig build -Denable-oracle=true upstream-oracle
```

`tools/check-build-policy`, the in-graph import test, and the standalone
`build_support/consumer` package enforce the initial dependency boundary. The
standalone fixture has its own test-only path manifest so it exercises the
same public package/module interface as an external project. CI also performs
an ordinary build with an empty package cache and a deliberately unusable
proxy, proving that the lazy oracle is not requested.
`tools/check-install-isolation` additionally inspects the installed static
library tree and rejects oracle-named files or pinned-upstream C++ symbols.

## Admission and removal

Every proposed dependency needs a Beads issue and all of:

1. immutable content hash and non-branch URL;
2. license and provenance review, including preserved notices;
3. evidence for every supported target;
4. a named owner and maintenance/removal plan;
5. justification against `std` or local code;
6. proof that its types do not escape through public APIs;
7. fetch-then-offline build evidence.

Test/tool-only dependencies must be lazy or isolated in a separate package and
must not be linked or imported by production artifacts. Remove a dependency
when its last justified consumer disappears; delete its manifest entry,
license/provenance record when legally permitted, build imports, cache-specific
tests, and admission Bead only after ordinary and offline package checks pass.
