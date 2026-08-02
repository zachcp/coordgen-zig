# Review probes (cgz-7v2 independent review, 2026-08-01)

Reusable harnesses behind the measurements in cgz-r05 / cgz-r13. Uncommitted; delete
freely — they are reference implementations for the Phase 0 conformance beads, not
project source.

- `probe2.cpp` — 7 realistic SMILES + a global `operator new` that hands out descending
  addresses, inverting every pointer comparison. Used for the optimization-level,
  fp-contract, and Guard Malloc checks.
- `probe3.cpp` — seeded random-graph corpus generator (xorshift64*) + the same allocator
  adversary + raw-f32 coordinate dump. Source of the architecture and allocator-order
  divergence numbers.

Both need upstream coordgenlibs at `d20e735d96480385b2e257522288004038a08cc9` in `cg/`:

```
git clone https://github.com/schrodinger/coordgenlibs.git cg && git -C cg checkout d20e735
zig c++ -std=c++17 -O2 -DIN_COORDGEN -I cg -o probe3_arm cg/*.cpp probe3.cpp
zig c++ -std=c++17 -O2 -target x86_64-macos -DIN_COORDGEN -I cg -o probe3_x64 cg/*.cpp probe3.cpp
./probe3_arm n 2000 > n.out; ./probe3_arm d 2000 > d.out; ./probe3_x64 n 2000 > x.out
diff n.out d.out | wc -l    #   112  -> 1/2000 molecules, max dev 738.9 units
diff n.out x.out | wc -l    # 64836  -> 1414/2000 molecules, max dev 2470.7 units
```

Bond length = 50 units.
