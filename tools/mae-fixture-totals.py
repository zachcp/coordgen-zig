#!/usr/bin/env python3
"""Derive the .mae fixture totals recorded in build.zig.

This is deliberately a second, independent implementation of the fixture
reading in src/conformance/mae.zig: it tokenises the file itself and resolves
the six upstream fields by name. Running it is how the expectations in
build.zig's `fixture_expectations` table are produced, so a build gate failure
means the two implementations disagree rather than that a digest drifted.

Usage: tools/mae-fixture-totals.py FIXTURE.mae...
"""

import re
import sys

TOKEN = re.compile(r'[^\s{}"#]+|[{}]')
INDEXED_BLOCK = re.compile(r"(m_atom|m_bond)\[(\d+)\]")


def tokenize(text):
    """Yield one token per value, dropping comments.

    No field this script reads is quoted, but a quoted value still occupies a
    column, so it is yielded whole - including the value in
    test/macrocycle.mae that spans twenty-one physical lines and the
    s_m_atom_name column inside templates.mae data rows.
    """
    index = 0
    while index < len(text):
        char = text[index]
        if char == '"':
            start = index
            index += 1
            while index < len(text) and text[index] != '"':
                index += 2 if text[index] == "\\" else 1
            index += 1
            yield text[start:index]
        elif char == "#":
            while index < len(text) and text[index] != "\n":
                index += 1
        elif char.isspace():
            index += 1
        else:
            match = TOKEN.match(text, index)
            yield match.group(0)
            index = match.end()


def totals(path):
    with open(path, encoding="utf-8") as handle:
        tokens = list(tokenize(handle.read()))

    result = dict(
        structures=0,
        atoms=0,
        bonds=0,
        atomic_number_sum=0,
        bond_order_sum=0,
        x_sum=0.0,
        y_sum=0.0,
    )

    index = 0
    while index < len(tokens):
        token = tokens[index]
        index += 1
        if token == "f_m_ct":
            result["structures"] += 1
            continue
        match = INDEXED_BLOCK.fullmatch(token)
        if match is None:
            continue

        kind, count = match.group(1), int(match.group(2))
        assert tokens[index] == "{", (path, token)
        index += 1
        properties = []
        while tokens[index] != ":::":
            properties.append(tokens[index])
            index += 1
        index += 1

        width = len(properties) + 1
        rows = [tokens[index + row * width:index + (row + 1) * width] for row in range(count)]
        index += count * width
        assert tokens[index] == ":::", (path, token, tokens[index])
        index += 1

        for ordinal, row in enumerate(rows, start=1):
            assert int(row[0]) == ordinal, (path, token, row[0])

        def column(name):
            return properties.index(name) + 1

        if kind == "m_atom":
            result["atoms"] += count
            for row in rows:
                result["atomic_number_sum"] += int(row[column("i_m_atomic_number")])
                result["x_sum"] += float(row[column("r_m_x_coord")])
                result["y_sum"] += float(row[column("r_m_y_coord")])
        else:
            result["bonds"] += count
            for row in rows:
                # Maestro atom references are 1-based.
                assert 1 <= int(row[column("i_m_from")]) <= result["atoms"], path
                assert 1 <= int(row[column("i_m_to")]) <= result["atoms"], path
                result["bond_order_sum"] += int(row[column("i_m_order")])

    return result


def main(paths):
    if not paths:
        print(__doc__.strip().splitlines()[-1], file=sys.stderr)
        return 2
    for path in paths:
        value = totals(path)
        print(
            "%s\tstructures=%d,atoms=%d,bonds=%d,atomic_number_sum=%d,"
            "bond_order_sum=%d,x_sum=%.6f,y_sum=%.6f"
            % (
                path,
                value["structures"],
                value["atoms"],
                value["bonds"],
                value["atomic_number_sum"],
                value["bond_order_sum"],
                value["x_sum"],
                value["y_sum"],
            )
        )
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
