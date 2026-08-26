# Validation audit

- Generated: 2026-08-25T18:57:42
- Commit: `0f715ef010ed9319f88eaab0519f52ce5ef57c5d` (WORKING TREE DIRTY)
- Toolchain: `0.17.0-dev.1516+8a4b5424d`
- Host: `arm64-darwin` — **one host only**; see the scope note at the end
- Reproduce: `tools/validation-audit`

Every row below was produced by running the command in it. A status is a measurement, not a description.

## Gates

| Step | Status | Reproduce |
|---|---|---|
| `install` | **NOT-A-GATE** | `installs artifacts; exercised by package-check and abi-check` |
| `uninstall` | **NOT-A-GATE** | `removes artifacts; not a gate` |
| `test` | **PASS** | `.tools/toolchains/aarch64-macos/0.17.0-dev.1516+8a4b5424d/zig build test` |
| `parity-ceiling-check` | **PASS** | `.tools/toolchains/aarch64-macos/0.17.0-dev.1516+8a4b5424d/zig build parity-ceiling-check` |
| `coverage-check` | **PASS** | `.tools/toolchains/aarch64-macos/0.17.0-dev.1516+8a4b5424d/zig build coverage-check` |
| `module-graph-check` | **PASS** | `.tools/toolchains/aarch64-macos/0.17.0-dev.1516+8a4b5424d/zig build module-graph-check` |
| `reachability-check` | **PASS** | `.tools/toolchains/aarch64-macos/0.17.0-dev.1516+8a4b5424d/zig build reachability-check` |
| `package-check` | **PASS** | `.tools/toolchains/aarch64-macos/0.17.0-dev.1516+8a4b5424d/zig build package-check` |
| `abi-check` | **PASS** | `.tools/toolchains/aarch64-macos/0.17.0-dev.1516+8a4b5424d/zig build abi-check` |
| `upstream-oracle` | **PASS** | `.tools/toolchains/aarch64-macos/0.17.0-dev.1516+8a4b5424d/zig build upstream-oracle -Denable-oracle=true` |
| `conformance` | **PASS** | `.tools/toolchains/aarch64-macos/0.17.0-dev.1516+8a4b5424d/zig build conformance -Denable-oracle=true` |
| `template-regeneration-check` | **PASS** | `.tools/toolchains/aarch64-macos/0.17.0-dev.1516+8a4b5424d/zig build template-regeneration-check -Denable-oracle=true` |
| `performance-baseline` | **PASS** | `.tools/toolchains/aarch64-macos/0.17.0-dev.1516+8a4b5424d/zig build performance-baseline -Denable-oracle=true` |
| `performance-check` | **PASS** | `.tools/toolchains/aarch64-macos/0.17.0-dev.1516+8a4b5424d/zig build performance-check -Denable-oracle=true` |
| `sanitizer-check` | **PASS** | `.tools/toolchains/aarch64-macos/0.17.0-dev.1516+8a4b5424d/zig build sanitizer-check -Denable-oracle=true -Dsanitize-oracle=true` |
| `examples` | **PASS** | `.tools/toolchains/aarch64-macos/0.17.0-dev.1516+8a4b5424d/zig build examples` |
| `fuzz` | **PASS** | `.tools/toolchains/aarch64-macos/0.17.0-dev.1516+8a4b5424d/zig build fuzz` |

## Reserved gates

A gate that fails by name rather than passing vacuously is what cgz-r21 requires. `reserved` means the work is not written yet and the message names its bead; `conditional` means the gate needs a flag or host capability and the message says which. Both are legitimate. Anything else is a gate that has been switched off.

| Kind | Owning bead | Message |
|---|---|---|
| conditional | - | sanitizer-check requires -Denable-oracle=true -Dsanitize-oracle=true |
| conditional | - | corpus-check needs to run {s}-{s} binaries for the architecture axis;  on Linux install qemu-user an |
| conditional | - | oracle steps require -Denable-oracle=true; ordinary builds never fetch conformance sources |

## Differential

Measured on the **drug_like** partition only, 7 members. Max coordinate deviation 0.092483 bond lengths.

| Observable | Matched |
|---|---|
| `clean_pose` | 7/7 |
| `coordinates` | 7/7 |
| `input_to_internal` | 7/7 |
| `internal_to_input` | 7/7 |
| `effective_bond_orders` | 7/7 |
| `bond_displays` | 7/7 |
| `atom_stereo` | 7/7 |

**Not measured against the oracle:** the `adversarial` partition (2000 members) is consumed by corpus-classify, not by native-oracle-diff, so no per-observable oracle comparison exists for it. No statement in this report covers inputs outside the `drug_like` partition.

## What this report does not establish

While any gate above is RESERVED or FAIL, the validation epic cannot close: an independent reviewer cannot reproduce a release gate that does not exist yet. This tool exits non-zero in that case by design.

**One host.** Every row was measured on `arm64-darwin`. The supported target matrix is wider, and nothing here says anything about the other legs of it - a gate that passes here may not have been run anywhere else. Completing that matrix is cgz-7v2.4.6.

**One optimize mode per run.** Gates were run at the build's default mode. Where a criterion names Debug AND ReleaseSafe, both have to be run, and this report covers whichever one this invocation used.

