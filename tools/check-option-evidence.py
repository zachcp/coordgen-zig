#!/usr/bin/env python3
"""Require executable evidence for every public generation option."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re
import sys
import tempfile


HEADER = ("option", "kind", "observable", "evidence", "decision")
KINDS = {"flip", "reject", "exempt"}


class EvidenceError(Exception):
    pass


@dataclass(frozen=True)
class Row:
    option: str
    kind: str
    observable: str
    evidence: str
    decision: str


def fail(message: str) -> None:
    raise EvidenceError(message)


def struct_fields(path: Path) -> tuple[str, ...]:
    source = path.read_text(encoding="utf-8")
    match = re.search(r"^pub const Options = (?:extern )?struct \{(.*?)^\};", source, re.MULTILINE | re.DOTALL)
    if match is None:
        fail(f"cannot find Options in {path}")
    return tuple(re.findall(r"^\s{4}([A-Za-z_][A-Za-z0-9_]*):", match.group(1), re.MULTILINE))


def read_rows(path: Path) -> list[Row]:
    lines = [line for line in path.read_text(encoding="utf-8").splitlines() if line and not line.startswith("#")]
    if not lines or tuple(lines[0].split("\t")) != HEADER:
        fail(f"manifest header must be {' TAB '.join(HEADER)}")
    rows: list[Row] = []
    seen: set[str] = set()
    for line_number, line in enumerate(lines[1:], start=2):
        fields = line.split("\t")
        if len(fields) != len(HEADER):
            fail(f"line {line_number}: expected {len(HEADER)} tab-separated fields")
        row = Row(*fields)
        if row.option in seen:
            fail(f"line {line_number}: duplicate option {row.option}")
        seen.add(row.option)
        if row.kind not in KINDS:
            fail(f"line {line_number}: kind must be flip, reject, or exempt")
        if not row.observable:
            fail(f"line {line_number}: option has no named observable or error")
        if row.kind == "exempt":
            if not row.decision.startswith("cgz-") or not row.evidence:
                fail(f"line {line_number}: exemption needs a reviewed Bead and reason")
        elif row.decision != "-":
            fail(f"line {line_number}: executable evidence must use '-' as its decision")
        rows.append(row)
    return rows


def evidence_exists(root: Path, row: Row) -> bool:
    if row.kind == "exempt":
        return True
    path_text, separator, needle = row.evidence.partition("::")
    path = root / path_text
    if not separator or not path.is_file() or not needle:
        return False
    source = path.read_text(encoding="utf-8")
    return needle in source and row.option in source


def validate(root: Path, manifest: Path) -> None:
    safe_fields = struct_fields(root / "src/api.zig")
    c_fields = struct_fields(root / "src/c_abi_types.zig")
    if tuple(field for field in c_fields if field != "reserved") != safe_fields:
        fail(f"C option mapping differs from safe API: safe={list(safe_fields)} C={list(c_fields)}")
    if c_fields.count("reserved") != 1:
        fail("C Options must contain exactly one non-option reserved slot")

    rows = read_rows(manifest)
    recorded = {row.option for row in rows}
    declared = set(safe_fields)
    missing = sorted(declared - recorded)
    stale = sorted(recorded - declared)
    if missing or stale:
        details = []
        if missing:
            details.append(f"missing={missing}")
        if stale:
            details.append(f"stale={stale}")
        fail("option evidence inventory differs: " + " ".join(details))
    for row in rows:
        if not evidence_exists(root, row):
            fail(f"{row.option}: evidence does not resolve or does not name the option: {row.evidence}")


def run_self_test(root: Path) -> None:
    with tempfile.TemporaryDirectory(prefix="cgz-option-evidence-") as directory:
        fixture = Path(directory)
        (fixture / "src").mkdir()
        (fixture / "tests").mkdir()
        (fixture / "src/api.zig").write_text(
            "pub const Options = struct {\n    alpha: bool = false,\n};\n", encoding="utf-8"
        )
        (fixture / "src/c_abi_types.zig").write_text(
            "pub const Options = extern struct {\n    alpha: u32 = 0,\n    reserved: u32 = 0,\n};\n",
            encoding="utf-8",
        )
        (fixture / "tests/evidence.zig").write_text(
            'test "alpha changes coordinates" { // alpha\n}\n', encoding="utf-8"
        )
        manifest = fixture / "manifest.tsv"
        manifest.write_text(
            "\t".join(HEADER)
            + "\nalpha\tflip\tcoordinates\ttests/evidence.zig::alpha changes coordinates\t-\n",
            encoding="utf-8",
        )
        validate(fixture, manifest)
        manifest.write_text("\t".join(HEADER) + "\n", encoding="utf-8")
        try:
            validate(fixture, manifest)
        except EvidenceError as error:
            if "missing=['alpha']" not in str(error):
                fail(f"self-test rejected the gap for the wrong reason: {error}")
        else:
            fail("self-test accepted a declared option with no evidence")


def main(argv: list[str]) -> int:
    root = Path(__file__).resolve().parent.parent
    try:
        if argv[1:] == ["--self-test"]:
            run_self_test(root)
            print("check-option-evidence: self-test rejected a planted uncovered option")
        elif argv[1:]:
            fail("usage: tools/check-option-evidence.py [--self-test]")
        else:
            validate(root, root / "conformance/option_evidence.tsv")
            print("check-option-evidence: every declared option has executable evidence or a reviewed exemption")
    except (EvidenceError, OSError) as error:
        print(f"check-option-evidence: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
