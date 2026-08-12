#!/usr/bin/env python3
"""Verify the public-test baseline and the pinned upstream harness defects."""

from __future__ import annotations

from collections import Counter
from pathlib import Path
import re
import sys


EXPECTED_CASES = {
    "test_coordgen.cpp": (
        "SampleTest",
        "TemplateTest",
        "ClearWedgesTest",
        "DisableMetalZOBs",
        "terminalMetalZOBs",
        "testMinimizedRingsShape",
        "testPolyominoCoordinatesOfSubstituent",
        "testPolyominoSameAs",
        "testGetDoubleBondConstraints",
        "testClockwiseOrderedSubstituents",
        "testClockwiseOrderedNaN",
        "testbicyclopentane",
        "testFusedRings",
        "testTemplates",
        "testRingComplex",
        "testCoordgenFragmenter",
    ),
    "test_smilesparser.cpp": ("Basics", "Rings", "Branching", "BondOrder"),
}

EXPECTED_ASSERTIONS = Counter(
    {
        "BOOST_CHECK": 8,
        "BOOST_CHECK_EQUAL": 1,
        "BOOST_CHECK_EQUAL_COLLECTIONS": 1,
        "BOOST_CHECK_MESSAGE": 2,
        "BOOST_REQUIRE": 40,
        "BOOST_REQUIRE_EQUAL": 19,
        "BOOST_TEST": 31,
    }
)

EXPECTED_FIXTURES = {
    "macrocycle.mae",
    "metalZobs.mae",
    "nonterminalMetalZobs.mae",
    "test.mae",
    "testChirality.mae",
    "test_mol.mae",
}


def fail(message: str) -> None:
    raise SystemExit(f"check-upstream-test-inventory: {message}")


def require(condition: bool, message: str) -> None:
    if not condition:
        fail(message)


def code_tokens(source: str) -> str:
    """Blank comments and literals so names inside them are not inventory."""
    ignored = re.compile(
        r"//[^\n]*|/\*.*?\*/|\"(?:\\.|[^\"\\])*\"|'(?:\\.|[^'\\])*'",
        re.DOTALL,
    )
    return ignored.sub(lambda match: re.sub(r"[^\n]", " ", match.group()), source)


def main() -> None:
    if len(sys.argv) != 2:
        fail("usage: check-upstream-test-inventory.py UPSTREAM_ROOT")

    root = Path(sys.argv[1])
    test_dir = root / "test"
    require(test_dir.is_dir(), f"test directory not found under {root}")

    assertions: Counter[str] = Counter()
    case_total = 0
    for filename, expected_cases in EXPECTED_CASES.items():
        source = (test_dir / filename).read_text(encoding="utf-8")
        tokens = code_tokens(source)
        cases = tuple(
            re.findall(r"\bBOOST_AUTO_TEST_CASE\s*\(\s*([A-Za-z_][A-Za-z0-9_]*)", tokens)
        )
        require(cases == expected_cases, f"{filename} test cases changed: {cases!r}")
        case_total += len(cases)
        assertions.update(
            re.findall(
                r"\b(BOOST_(?:CHECK|REQUIRE|TEST)(?:_[A-Z]+)*)\s*\(", tokens
            )
        )

    require(case_total == 20, f"expected 20 public test cases, found {case_total}")
    require(
        assertions == EXPECTED_ASSERTIONS,
        f"assertion macro inventory changed: {dict(sorted(assertions.items()))}",
    )
    require(sum(assertions.values()) == 102, "expected 102 assertion macro sites")

    fixtures = {path.name for path in test_dir.glob("*.mae")}
    require(fixtures == EXPECTED_FIXTURES, f"fixture inventory changed: {sorted(fixtures)}")

    test_cmake = (test_dir / "CMakeLists.txt").read_text(encoding="utf-8")
    require(
        re.search(
            r"add_test\s*\(\s*NAME\s+test_smilesparser\s+COMMAND\s+"
            r"\$\{CMAKE_BINARY_DIR\}/test/test_coordgen\b",
            test_cmake,
            re.DOTALL,
        )
        is not None,
        "expected pinned CTest miswiring is absent; re-audit the correction policy",
    )
    require(
        test_cmake.lstrip().startswith(
            "# Use the example executable as a test\nadd_test(example "
        ),
        "expected unguarded example registration is absent; re-audit the correction policy",
    )

    example = (root / "example_dir" / "example.cpp").read_text(encoding="utf-8")
    require(
        not re.search(r"\b(?:BOOST_(?:CHECK|REQUIRE|TEST)|assert)\s*\(", example),
        "upstream example gained an assertion; re-audit the local output contract",
    )

    root_cmake = (root / "CMakeLists.txt").read_text(encoding="utf-8")
    require(
        'set(MAEPARSER_VERSION "master"' in root_cmake,
        "moving maeparser dependency pin changed; re-audit dependency policy",
    )
    require(
        "MEMORYCHECK_COMMAND_OPTIONS" in root_cmake,
        "upstream memory-check configuration changed; re-audit the local memory gates",
    )
    pipeline_text = "\n".join(
        path.read_text(encoding="utf-8")
        for path in sorted((root / ".azure-pipelines").glob("*.yml"))
    ) + (root / "azure-pipelines.yml").read_text(encoding="utf-8")
    require(
        re.search(r"ctest\s+(?:-T\s+MemCheck|--test-action\s+MemCheck)", pipeline_text)
        is None,
        "upstream CI gained a CTest memory-check action; re-audit the local policy",
    )

    smiles_header = (test_dir / "coordgenBasicSMILES.h").read_text(encoding="utf-8")
    require(
        'sketcherMinimizerAtom.h' not in smiles_header,
        "SMILES helper became self-contained; re-audit the include-order correction",
    )

    details = ", ".join(
        f"{name.removeprefix('BOOST_')}={count}"
        for name, count in sorted(assertions.items())
    )
    print(
        "verified pinned upstream public-test inventory: "
        f"20 cases, 102 assertion sites ({details}), 6 fixtures, 1 example; "
        "all six documented harness defects remain present upstream and corrected locally"
    )


if __name__ == "__main__":
    main()
