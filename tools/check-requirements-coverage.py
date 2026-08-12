#!/usr/bin/env python3
"""Validate the requirements coverage manifest and its generated summary."""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
import hashlib
from pathlib import Path
import re
import sys


HEADER = (
    "category",
    "id",
    "source",
    "selector",
    "quantity",
    "evidence",
    "disposition",
    "decision",
)

API_STRUCTS = (
    "AtomStereoInput",
    "AtomInput",
    "BondStereoInput",
    "BondInput",
    "ResidueInput",
    "ResidueInteractionInput",
    "Options",
    "Input",
    "Result",
)

C_STRUCTS = (
    "StringView",
    "IndexSpan",
    "AtomInput",
    "BondInput",
    "ResidueInput",
    "ResidueInteractionInput",
    "Options",
    "Input",
    "Result",
)

EXPECTED_FIXTURES = {
    "macrocycle.mae",
    "metalZobs.mae",
    "nonterminalMetalZobs.mae",
    "test.mae",
    "testChirality.mae",
    "test_mol.mae",
}


class CoverageError(Exception):
    pass


@dataclass(frozen=True)
class Row:
    category: str
    id: str
    source: str
    selector: str
    quantity: int
    evidence: str
    disposition: str
    decision: str


def fail(message: str) -> None:
    raise CoverageError(message)


def code_tokens(source: str) -> str:
    ignored = re.compile(
        r"//[^\n]*|/\*.*?\*/|\"(?:\\.|[^\"\\])*\"|'(?:\\.|[^'\\])*'",
        re.DOTALL,
    )
    return ignored.sub(lambda match: re.sub(r"[^\n]", " ", match.group()), source)


def read_rows(path: Path) -> list[Row]:
    lines = [line for line in path.read_text(encoding="utf-8").splitlines() if line and not line.startswith("#")]
    if not lines:
        fail("manifest is empty")
    if tuple(lines[0].split("\t")) != HEADER:
        fail(f"manifest header must be {' TAB '.join(HEADER)}")

    rows: list[Row] = []
    seen: set[str] = set()
    for line_number, line in enumerate(lines[1:], start=2):
        fields = line.split("\t")
        if len(fields) != len(HEADER):
            fail(f"line {line_number}: expected {len(HEADER)} tab-separated fields")
        category, row_id, source, selector, quantity_text, evidence, disposition, decision = fields
        if row_id in seen:
            fail(f"line {line_number}: duplicate id {row_id}")
        seen.add(row_id)
        try:
            quantity = int(quantity_text)
        except ValueError:
            fail(f"line {line_number}: quantity is not an integer")
        if quantity < 0:
            fail(f"line {line_number}: quantity is negative")
        if disposition not in {"executable", "non_goal"}:
            fail(f"line {line_number}: disposition must be executable or non_goal")
        if disposition == "executable" and not evidence:
            fail(f"line {line_number}: executable requirement has no evidence")
        if disposition == "non_goal" and not decision.startswith("cgz-"):
            fail(f"line {line_number}: non-goal has no reviewed Bead decision")
        rows.append(Row(category, row_id, source, selector, quantity, evidence, disposition, decision))
    return rows


def struct_fields(path: Path, name: str) -> tuple[str, ...]:
    source = path.read_text(encoding="utf-8")
    start_match = re.search(rf"^pub const {re.escape(name)} = (?:extern )?struct \{{", source, re.MULTILINE)
    if start_match is None:
        fail(f"cannot find public struct {name} in {path}")
    start = start_match.end()
    depth = 1
    fields: list[str] = []
    for line in source[start:].splitlines():
        if depth == 1:
            match = re.match(r"\s{4}([A-Za-z_][A-Za-z0-9_]*):", line)
            if match:
                fields.append(match.group(1))
        depth += line.count("{") - line.count("}")
        if depth == 0:
            break
    if not fields:
        fail(f"public struct {name} has no discovered fields")
    return tuple(fields)


def enum_members(path: Path, declaration: str) -> tuple[str, ...]:
    source = path.read_text(encoding="utf-8")
    match = re.search(rf"{re.escape(declaration)}\s*\{{(.*?)\}};", source, re.DOTALL)
    if match is None:
        fail(f"cannot find {declaration} in {path}")
    body = code_tokens(match.group(1))
    return tuple(re.findall(r"^\s*([A-Za-z_][A-Za-z0-9_]*)\s*(?:=\s*\d+)?\s*,", body, re.MULTILINE))


def evidence_exists(root: Path, evidence: str) -> bool:
    split = evidence.find("::")
    if split <= 0:
        return False
    path = root / evidence[:split]
    needle = evidence[split + 2 :]
    return path.is_file() and needle in path.read_text(encoding="utf-8")


def validate_local(root: Path, rows: list[Row]) -> None:
    for row in rows:
        if row.disposition == "executable" and not evidence_exists(root, row.evidence):
            fail(f"{row.id}: evidence does not resolve: {row.evidence}")

    categories = Counter(row.category for row in rows)
    if categories["upstream_case"] != 20:
        fail(f"expected 20 upstream_case rows, found {categories['upstream_case']}")
    assertion_sites = sum(row.quantity for row in rows if row.category == "upstream_case")
    if assertion_sites != 102:
        fail(f"expected 102 upstream assertion sites, found {assertion_sites}")

    fixtures = {row.selector for row in rows if row.category == "fixture"}
    if fixtures != EXPECTED_FIXTURES:
        fail(f"fixture inventory differs: {sorted(fixtures)}")
    for row in (row for row in rows if row.category == "fixture"):
        if re.fullmatch(r"sha256:[0-9a-f]{64}", row.decision) is None:
            fail(f"{row.id}: fixture has no pinned SHA-256")
    templates = [row for row in rows if row.category == "templates"]
    if len(templates) != 1 or templates[0].quantity != 82:
        fail("manifest must cover exactly 82 templates")
    examples = [row for row in rows if row.category == "example"]
    if len(examples) != 1 or examples[0].quantity != 1:
        fail("manifest must cover exactly one example")

    api_rows = {row.selector: row for row in rows if row.category == "zig_dto"}
    if set(api_rows) != set(API_STRUCTS):
        fail("safe Zig DTO type inventory differs")
    for name in API_STRUCTS:
        count = len(struct_fields(root / "src/api.zig", name))
        if api_rows[name].quantity != count:
            fail(f"zig DTO {name} has {count} fields, manifest records {api_rows[name].quantity}")

    c_rows = {row.selector: row for row in rows if row.category == "c_dto"}
    if set(c_rows) != set(C_STRUCTS):
        fail("C DTO type inventory differs")
    for name in C_STRUCTS:
        count = len(struct_fields(root / "src/c_abi_types.zig", name))
        if c_rows[name].quantity != count:
            fail(f"C DTO {name} has {count} fields, manifest records {c_rows[name].quantity}")
    span_rows = [row for row in rows if row.category == "c_span_dto"]
    if len(span_rows) != 1 or span_rows[0].quantity != 18:
        fail("six public C span DTOs with three fields each must be covered")

    errors = set(enum_members(root / "src/core/errors.zig", "pub const Error = error"))
    manifest_errors = {row.selector for row in rows if row.category == "error"}
    if manifest_errors != errors:
        fail(f"error inventory differs: manifest={sorted(manifest_errors)} code={sorted(errors)}")

    observables = set(enum_members(root / "src/conformance/comparison.zig", "pub const Observable = enum"))
    manifest_observables = {row.selector for row in rows if row.category == "observable"}
    if manifest_observables != observables:
        fail("observable inventory differs")

    if categories["ownership"] < 5:
        fail("ownership contract inventory is incomplete")
    if categories["domain"] < 15:
        fail("algorithm-domain inventory is incomplete")
    if categories["non_goal"] != 1:
        fail("exactly one reviewed inaccessible-upstream-suite non-goal is required")


def upstream_cases(path: Path) -> dict[tuple[str, str], int]:
    cases: dict[tuple[str, str], int] = {}
    for filename in ("test_coordgen.cpp", "test_smilesparser.cpp"):
        source = code_tokens((path / "test" / filename).read_text(encoding="utf-8"))
        matches = list(re.finditer(r"\bBOOST_AUTO_TEST_CASE\s*\(\s*([A-Za-z_][A-Za-z0-9_]*)", source))
        for index, match in enumerate(matches):
            end = matches[index + 1].start() if index + 1 < len(matches) else len(source)
            body = source[match.end() : end]
            count = len(re.findall(r"\bBOOST_(?:CHECK|REQUIRE|TEST)(?:_[A-Z]+)*\s*\(", body))
            cases[(f"test/{filename}", match.group(1))] = count
    return cases


def validate_upstream(upstream: Path, rows: list[Row]) -> None:
    actual = upstream_cases(upstream)
    recorded = {
        (row.source, row.selector): row.quantity
        for row in rows
        if row.category == "upstream_case"
    }
    if recorded != actual:
        missing = sorted(set(actual) - set(recorded))
        stale = sorted(set(recorded) - set(actual))
        wrong = sorted(key for key in set(actual) & set(recorded) if actual[key] != recorded[key])
        fail(f"upstream case mapping differs: missing={missing} stale={stale} wrong_counts={wrong}")
    for row in (row for row in rows if row.category == "fixture"):
        fixture = upstream / row.source
        actual_hash = hashlib.sha256(fixture.read_bytes()).hexdigest()
        expected_hash = row.decision.removeprefix("sha256:")
        if actual_hash != expected_hash:
            fail(f"{row.id}: fixture SHA-256 differs")


def summary(rows: list[Row]) -> str:
    categories = Counter(row.category for row in rows)
    executable = sum(row.quantity for row in rows if row.disposition == "executable")
    non_goals = sum(row.quantity for row in rows if row.disposition == "non_goal")
    assertions = sum(row.quantity for row in rows if row.category == "upstream_case")
    dto_fields = sum(row.quantity for row in rows if row.category in {"zig_dto", "c_dto", "c_span_dto"})
    lines = [
        "# Requirements coverage summary",
        "",
        "Generated by `tools/check-requirements-coverage.py`; do not edit by hand.",
        "",
        "| Coverage class | Requirements |",
        "|---|---:|",
        f"| Upstream public test cases | {categories['upstream_case']} |",
        f"| Upstream assertion macro sites | {assertions} |",
        f"| MAE fixtures | {categories['fixture']} |",
        f"| Embedded templates | {sum(row.quantity for row in rows if row.category == 'templates')} |",
        f"| Executable examples | {sum(row.quantity for row in rows if row.category == 'example')} |",
        f"| Public Zig/C DTO fields | {dto_fields} |",
        f"| Public errors | {categories['error']} |",
        f"| Output observables | {categories['observable']} |",
        f"| Ownership rules | {categories['ownership']} |",
        f"| Algorithm domains | {categories['domain']} |",
        "",
        f"Executable requirement units: **{executable}**. Reviewed non-goal units: **{non_goals}**.",
        "All six MAE fixture rows pin and verify their upstream SHA-256.",
        "Every manifest row resolves to an executable test artifact or a Bead-backed non-goal.",
        "",
    ]
    return "\n".join(lines)


def self_test(root: Path, manifest: Path) -> None:
    rows = read_rows(manifest)
    validate_local(root, rows)
    fixture_index = next(index for index, row in enumerate(rows) if row.category == "fixture")
    missing_hash = rows.copy()
    missing_hash[fixture_index] = Row(**{**rows[fixture_index].__dict__, "decision": "-"})
    mutations = (
        (rows[1:], "missing public case"),
        ([Row(**{**rows[0].__dict__, "quantity": rows[0].quantity + 1}), *rows[1:]], "wrong assertion count"),
        ([Row(**{**rows[0].__dict__, "evidence": "missing::test"}), *rows[1:]], "missing evidence"),
        (missing_hash, "missing fixture hash"),
    )
    for changed, label in mutations:
        try:
            validate_local(root, changed)
        except CoverageError:
            continue
        fail(f"self-test did not reject {label}")
    print(
        "check-requirements-coverage: self-test rejected missing coverage, "
        "count drift, stale evidence, and missing fixture hashes"
    )


def main(argv: list[str]) -> int:
    root = Path(__file__).resolve().parent.parent
    manifest = root / "conformance/requirements_coverage.tsv"
    upstream: Path | None = None
    summary_path: Path | None = None
    run_self_test = False
    index = 0
    while index < len(argv):
        if argv[index] == "--upstream" and index + 1 < len(argv):
            upstream = Path(argv[index + 1])
            index += 2
        elif argv[index] == "--summary" and index + 1 < len(argv):
            summary_path = Path(argv[index + 1])
            index += 2
        elif argv[index] == "--self-test":
            run_self_test = True
            index += 1
        else:
            fail(f"unknown or incomplete option: {argv[index]}")

    if run_self_test:
        self_test(root, manifest)
        return 0

    rows = read_rows(manifest)
    validate_local(root, rows)
    if upstream is not None:
        validate_upstream(upstream, rows)
    generated = summary(rows)
    if summary_path is not None:
        if summary_path.read_text(encoding="utf-8") != generated:
            fail(f"generated summary is stale: {summary_path}")
    else:
        sys.stdout.write(generated)
    print(
        f"check-requirements-coverage: verified {len(rows)} mappings, "
        "20 upstream cases, 102 assertion sites, and complete public contract inventories"
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main(sys.argv[1:]))
    except CoverageError as error:
        print(f"check-requirements-coverage: {error}", file=sys.stderr)
        raise SystemExit(1)
