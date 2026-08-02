# Toolchain policy

`toolchain.lock` is the pre-execution trust root. It is a small tab-separated,
host-readable file rather than ZON because the compiler must be authenticated
before Zig is allowed to parse `build.zig.zon`. Each supported host records the
exact distribution archive name, byte length, SHA-256, immutable upstream URL,
and owner-controlled mirror URL.

`tools/bootstrap-zig` detects the host, downloads the selected distribution,
checks its byte length and SHA-256 before extraction, and then checks the
extracted compiler's exact `zig version`. CI obtains Zig only through this
script and invokes no preinstalled `zig` first.

The supported build contract intentionally targets
`0.17.0-dev.1516+8a4b5424d` and its new Maker/package behavior. Whole-project
0.16 compatibility is not a gate. Prefer the tagged 0.17 release once it
exists, using the bump procedure below.

## Mirror requirement

The mirror column is currently `UNCONFIGURED` because this repository has no
remote and no owner-controlled artifact endpoint was supplied. The bootstrap
therefore warns and uses the exact checksummed Zig CDN URL. Set
`CGZ_REQUIRE_TOOLCHAIN_MIRROR=1` to make this missing external input fatal.
Replace each marker with an immutable URL under project control before closing
the toolchain provenance issue; the existing SHA-256 remains authoritative.

## Bump procedure

1. Create or identify the Beads toolchain-bump issue.
2. Update the exact version and every supported-host archive URL, byte length,
   and SHA-256 in `toolchain.lock`. Mirror those same bytes under project
   control and update the mirror URLs.
3. Update only `minimum_zig_version` in `build.zig.zon`; compiler distribution
   metadata never belongs in ZON.
4. Bootstrap the new compiler and run all build, test, package, ABI, template,
   fuzz, and conformance gates on each supported host.
5. Regenerate architecture/toolchain/optimization-keyed oracle baselines, diff
   them against the previous pin, and attach the complete diff to the bump
   Bead. A compiler bump does not change the semantic CoordGen pin.
6. Run `bd lint` and record the commands/results on the bump Bead before it is
   closed.
