#!/usr/bin/env python3
"""Independently extract and normalize upstream CoordgenTemplates.cpp."""
import math
import re
import struct
import sys


def f32(value):
    return struct.unpack("<f", struct.pack("<f", value))[0]


def bits(value):
    return struct.unpack("<I", struct.pack("<f", value))[0]


def normalize(atoms, bonds):
    dds, counts = [], []
    for a, b, _ in bonds:
        dx = f32(atoms[a][1] - atoms[b][1])
        dy = f32(atoms[a][2] - atoms[b][2])
        dd = f32(f32(dx * dx) + f32(dy * dy))
        for i, representative in enumerate(dds):
            if dd * .9 < representative and dd * 1.1 > representative:
                counts[i] += 1
                break
        else:
            dds.append(dd)
            counts.append(1)
    if not dds:
        return atoms
    selected = max(range(len(counts)), key=lambda i: counts[i])
    scale = f32(math.sqrt(dds[selected]))
    return [(z, f32(x / scale), f32(y / scale)) for z, x, y in atoms]


def main():
    if len(sys.argv) not in (3, 4) or (len(sys.argv) == 4 and sys.argv[3] != "--raw"):
        raise SystemExit("usage: extract-template-reference.py CoordgenTemplates.cpp OUTPUT [--raw]")
    raw = len(sys.argv) == 4
    text = open(sys.argv[1], encoding="utf-8").read()
    blocks = re.findall(r"std::array<std::tuple<int, float, float>.*?atoms = \{\{(.*?)\}\};\s*std::array<std::array<int, 3>.*?bonds = \{\{(.*?)\}\};", text, re.S)
    templates = []
    for atom_text, bond_text in blocks:
        atoms = [(int(z), f32(float(x)), f32(float(y))) for z, x, y in re.findall(r"tuple<int, float, float>\((\d+),\s*([-+\d.eE]+)f,\s*([-+\d.eE]+)f\)", atom_text)]
        bonds = [tuple(map(int, row)) for row in re.findall(r"\{\s*(\d+),\s*(\d+),\s*(\d+)\s*\}", bond_text)]
        templates.append((atoms if raw else normalize(atoms, bonds), bonds))
    if (len(templates), sum(len(a) for a, _ in templates), sum(len(b) for _, b in templates)) != (82, 1704, 1963):
        raise SystemExit("unexpected C++ template totals")
    out = open(sys.argv[2], "w", encoding="ascii", newline="\n")
    out.write("// Generated template reference: do not edit.\n")
    out.write("pub const Atom = struct { atomic_number: u8, x_bits: u32, y_bits: u32 };\n")
    out.write("pub const Bond = struct { from: u8, to: u8, order: u8 };\n")
    out.write("pub const Template = struct { atom_start: u16, atom_len: u8, bond_start: u16, bond_len: u8 };\n")
    out.write("pub const atoms = [_]Atom{\n")
    for atoms, _ in templates:
        for z, x, y in atoms:
            out.write(f"    .{{ .atomic_number = {z}, .x_bits = 0x{bits(x):08x}, .y_bits = 0x{bits(y):08x} }},\n")
    out.write("};\npub const bonds = [_]Bond{\n")
    for _, bonds in templates:
        for a, b, order in bonds:
            out.write(f"    .{{ .from = {a}, .to = {b}, .order = {order} }},\n")
    out.write("};\npub const templates = [_]Template{\n")
    ai = bi = 0
    for atoms, bonds in templates:
        out.write(f"    .{{ .atom_start = {ai}, .atom_len = {len(atoms)}, .bond_start = {bi}, .bond_len = {len(bonds)} }},\n")
        ai += len(atoms); bi += len(bonds)
    out.write("};\n")
    out.close()


if __name__ == "__main__":
    main()
