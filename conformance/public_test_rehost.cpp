/*
 * Boost-free rehosting of the pinned public tests that do not require MAE
 * fixtures. Test names and assertions retain their upstream provenance in
 * conformance/requirements_coverage.tsv. This is conformance-only code and is
 * never installed.
 */

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "CoordgenFragmenter.h"
#include "CoordgenMacrocycleBuilder.h"
#include "sketcherMinimizer.h"
#include "sketcherMinimizerAtom.h"
#include "sketcherMinimizerBond.h"
#include "sketcherMinimizerMolecule.h"
#include "sketcherMinimizerMaths.h"
#include "test/coordgenBasicSMILES.h"

namespace
{

using schrodinger::approxSmilesParse;

void check(bool condition, const char* message)
{
    if (!condition) {
        throw std::runtime_error(message);
    }
}

std::map<sketcherMinimizerAtom*, int>
reportingIndices(sketcherMinimizerMolecule& molecule)
{
    std::map<sketcherMinimizerAtom*, int> indices;
    int index = 0;
    for (auto* atom : molecule.getAtoms()) {
        indices.emplace(atom, ++index);
    }
    return indices;
}

bool bondsNearIdeal(sketcherMinimizerMolecule& molecule,
                    std::map<sketcherMinimizerAtom*, int>& indices,
                    const std::set<std::pair<int, int>>& skip = {})
{
    const float target = BONDLENGTH * BONDLENGTH;
    const float tolerance = target * 0.1f;
    for (auto* bond : molecule.getBonds()) {
        const std::pair<int, int> pair(indices[bond->getStartAtom()],
                                       indices[bond->getEndAtom()]);
        if (skip.count(pair) != 0) {
            continue;
        }
        const float distance = sketcherMinimizerMaths::squaredDistance(
            bond->getStartAtom()->getCoordinates(),
            bond->getEndAtom()->getCoordinates());
        if (distance < target - tolerance || distance > target + tolerance) {
            return false;
        }
    }
    return true;
}

bool noCrossingBonds(sketcherMinimizerMolecule& molecule)
{
    for (auto* first : molecule.getBonds()) {
        for (auto* second : molecule.getBonds()) {
            if (first == second ||
                first->getStartAtom() == second->getStartAtom() ||
                first->getStartAtom() == second->getEndAtom() ||
                first->getEndAtom() == second->getStartAtom() ||
                first->getEndAtom() == second->getEndAtom()) {
                continue;
            }
            if (sketcherMinimizerMaths::intersectionOfSegments(
                    first->getStartAtom()->getCoordinates(),
                    first->getEndAtom()->getCoordinates(),
                    second->getStartAtom()->getCoordinates(),
                    second->getEndAtom()->getCoordinates())) {
                return false;
            }
        }
    }
    return true;
}

void destroyParsedMolecule(sketcherMinimizerMolecule* molecule)
{
    for (auto* bond : molecule->getBonds()) {
        delete bond;
    }
    for (auto* atom : molecule->getAtoms()) {
        delete atom;
    }
    delete molecule;
}

void Basics()
{
    auto* molecule = approxSmilesParse("CCCC");
    check(molecule->getAtoms().size() == 4, "Basics atom count");
    destroyParsedMolecule(molecule);

    molecule = approxSmilesParse("CNO");
    const auto& atoms = molecule->getAtoms();
    auto* carbon = atoms[0];
    auto* nitrogen = atoms[1];
    auto* oxygen = atoms[2];
    check(carbon->getAtomicNumber() == 6, "Basics carbon");
    check(nitrogen->getAtomicNumber() == 7, "Basics nitrogen");
    check(oxygen->getAtomicNumber() == 8, "Basics oxygen");
    check(carbon->isNeighborOf(nitrogen), "Basics C-N neighbor");
    check(nitrogen->isNeighborOf(oxygen), "Basics N-O neighbor");
    check(!carbon->isNeighborOf(oxygen), "Basics C-O non-neighbor");
    destroyParsedMolecule(molecule);
}

void Rings()
{
    auto* molecule = approxSmilesParse("C1CCC1C");
    check(molecule->getAtoms().size() == 5, "Rings atom count");
    check(molecule->getBonds().size() == 5, "Rings bond count");
    check(molecule->getAtoms()[0]->isNeighborOf(molecule->getAtoms()[3]),
          "Rings closure");
    destroyParsedMolecule(molecule);
}

void Branching()
{
    auto* molecule = approxSmilesParse("CC(C)(C)C");
    check(molecule->getAtoms()[1]->getBonds().size() == 4,
          "Branching quaternary center");
    destroyParsedMolecule(molecule);

    molecule = approxSmilesParse("CC(C)(CC)C");
    check(molecule->getAtoms()[3]->getBonds().size() == 2,
          "Branching chain atom");
    destroyParsedMolecule(molecule);
}

void BondOrder()
{
    auto* molecule = approxSmilesParse("C=C");
    check(molecule->getBonds()[0]->getBondOrder() == 2,
          "BondOrder double bond");
    destroyParsedMolecule(molecule);
}

void testPolyominoCoordinatesOfSubstituent()
{
    Polyomino polyomino;
    polyomino.addHex(hexCoords(0, 0));
    check(polyomino.coordinatesOfSubstituent(vertexCoords(1, 0, 0)) ==
              vertexCoords(1, -1, -1),
          "polyomino one hex substituent");
    polyomino.addHex(hexCoords(1, 0));
    check(polyomino.coordinatesOfSubstituent(vertexCoords(1, 0, 0)) ==
              vertexCoords(0, 0, -1),
          "polyomino two hex substituent");
}

Polyomino polyomino(std::initializer_list<hexCoords> coordinates)
{
    Polyomino result;
    for (const auto coordinate : coordinates) {
        result.addHex(coordinate);
    }
    return result;
}

void testPolyominoSameAs()
{
    auto first = polyomino({{0, 0}, {1, 0}, {2, 0}, {0, 1}});
    check(first.isTheSameAs(first), "polyomino identity");

    auto reordered = polyomino({{2, 0}, {0, 0}, {0, 1}, {1, 0}});
    check(first.isTheSameAs(reordered), "polyomino order independence");
    check(reordered.isTheSameAs(first), "polyomino reverse order independence");

    auto translated = polyomino({{4, 2}, {5, 2}, {6, 2}, {4, 3}});
    check(first.isTheSameAs(translated), "polyomino translation");
    check(translated.isTheSameAs(first), "polyomino reverse translation");

    auto rotated = polyomino({{0, 0}, {-1, 0}, {-2, 0}, {0, -1}});
    check(first.isTheSameAs(rotated), "polyomino rotation");
    check(rotated.isTheSameAs(first), "polyomino reverse rotation");

    auto reflected = polyomino({{0, 0}, {0, 1}, {0, 2}, {1, 0}});
    check(!first.isTheSameAs(reflected), "polyomino rejects reflection");
    check(!reflected.isTheSameAs(first), "polyomino reverse reflection");

    auto shorter = polyomino({{1, 0}, {2, 0}, {0, 1}});
    check(!first.isTheSameAs(shorter), "polyomino rejects shorter set");
    check(!shorter.isTheSameAs(first), "polyomino reverse shorter set");
}

void testClockwiseOrderedSubstituents()
{
    auto* molecule = approxSmilesParse("CN(C)C");
    sketcherMinimizer minimizer;
    minimizer.initialize(molecule);
    minimizer.runGenerateCoordinates();

    const auto& atoms = minimizer.getAtoms();
    auto* center = atoms.at(0);
    auto* first = atoms.at(1);
    auto* second = atoms.at(2);
    auto* third = atoms.at(3);
    check(center->getAtomicNumber() == 7, "clockwise center nitrogen");

    const sketcherMinimizerPointF above(0, 50);
    const sketcherMinimizerPointF left(-50, 0);
    const sketcherMinimizerPointF right(50, 0);
    center->coordinates = {0, 0};
    first->coordinates = above;
    second->coordinates = left;
    third->coordinates = right;
    auto ordered = center->clockwiseOrderedNeighbors();
    check(ordered[0] == first && ordered[1] == second && ordered[2] == third,
          "clockwise first ordering");
    check((ordered[0]->coordinates - above).length() == 0 &&
              (ordered[1]->coordinates - left).length() == 0 &&
              (ordered[2]->coordinates - right).length() == 0,
          "clockwise first coordinates");

    third->coordinates = left;
    second->coordinates = right;
    ordered = center->clockwiseOrderedNeighbors();
    check(ordered[0] == first && ordered[1] == third && ordered[2] == second,
          "clockwise second ordering");
    check((ordered[0]->coordinates - above).length() == 0 &&
              (ordered[1]->coordinates - left).length() == 0 &&
              (ordered[2]->coordinates - right).length() == 0,
          "clockwise second coordinates");
}

void testClockwiseOrderedNaN()
{
    auto* molecule = approxSmilesParse("CN(C)C");
    molecule->getAtoms().at(1)->coordinates = {
        std::nanf("coordgen"), std::nanf("coordgen")};
    (void)molecule->getAtoms().at(0)->clockwiseOrderedNeighbors();
    destroyParsedMolecule(molecule);
}

void testbicyclopentane()
{
    auto* molecule = approxSmilesParse("C1C2CC1C2");
    sketcherMinimizer minimizer;
    minimizer.initialize(molecule);
    minimizer.runGenerateCoordinates();
    const auto& atoms = minimizer.getAtoms();
    auto* first = atoms.at(1);
    auto* second = atoms.at(2);
    auto* third = atoms.at(3);
    check(first->neighbors.size() == 2, "bicyclopentane first degree");
    check(second->neighbors.size() == 2, "bicyclopentane second degree");
    check(third->neighbors.size() == 2, "bicyclopentane third degree");
    check((first->getCoordinates() - second->getCoordinates()).length() > 15,
          "bicyclopentane first distance");
    check((first->getCoordinates() - third->getCoordinates()).length() > 15,
          "bicyclopentane second distance");
    check((second->getCoordinates() - third->getCoordinates()).length() > 15,
          "bicyclopentane third distance");
}

void testFusedRings()
{
    const std::vector<std::string> smiles{
        "C1CCC23CCCCC2CC3C1",
        "C1=CC2C3CC4C(CC3NC3CCCC(C32)N1)NC1CCCC2C1C4CCN2",
    };
    for (const auto& value : smiles) {
        auto* molecule = approxSmilesParse(value);
        sketcherMinimizer minimizer;
        minimizer.initialize(molecule);
        minimizer.runGenerateCoordinates();
        auto indices = reportingIndices(*molecule);
        check(bondsNearIdeal(*molecule, indices), "fused ring bond lengths");
        check(noCrossingBonds(*molecule), "fused ring crossings");
    }
}

void testTemplates()
{
    auto* molecule = approxSmilesParse("C12CC3CC(CC3C1)C2");
    sketcherMinimizer minimizer;
    minimizer.initialize(molecule);
    minimizer.runGenerateCoordinates();
    auto indices = reportingIndices(*molecule);
    const std::set<std::pair<int, int>> skip{{2, 5}, {6, 2}};
    check(bondsNearIdeal(*molecule, indices, skip), "template bond lengths");
}

void testRingComplex()
{
    auto* molecule = approxSmilesParse("CC1CC2CCCC(C3CC4CCC(C3)C4C)C2O1");
    sketcherMinimizer minimizer;
    minimizer.initialize(molecule);
    minimizer.runGenerateCoordinates();
    check(noCrossingBonds(*molecule), "ring complex crossings");
}

void testCoordgenFragmenter()
{
    auto* molecule = approxSmilesParse("CCCCC1CC(CCC)C1");
    const auto atoms = molecule->getAtoms();
    for (int index = 3; index <= 8; ++index) {
        atoms[index - 1]->constrained = true;
    }
    atoms[10]->constrained = true;
    const auto atom_indices = reportingIndices(*molecule);

    sketcherMinimizer minimizer;
    minimizer.initialize(molecule);
    CoordgenFragmenter::splitIntoFragments(molecule);
    const std::vector<std::set<int>> expected{
        {1, 2}, {3}, {4}, {5, 6, 7, 11}, {8}, {9, 10},
    };
    std::vector<std::set<int>> actual;
    for (auto* fragment : molecule->_fragments) {
        std::set<int> members;
        for (auto* atom : fragment->getAtoms()) {
            members.insert(atom_indices.at(atom));
        }
        actual.push_back(std::move(members));
    }
    std::sort(actual.begin(), actual.end());
    check(actual == expected, "fragment membership");
    check(!atoms[0]->fragment->constrained &&
              !atoms[0]->fragment->constrainedFlip,
          "fragment 1 flags");
    check(atoms[2]->fragment->constrained &&
              !atoms[2]->fragment->constrainedFlip,
          "fragment 2 flags");
    check(atoms[3]->fragment->constrained &&
              atoms[3]->fragment->constrainedFlip,
          "fragment 3 flags");
    check(atoms[4]->fragment->constrained &&
              atoms[4]->fragment->constrainedFlip,
          "fragment 4 flags");
    check(atoms[7]->fragment->constrained &&
              !atoms[7]->fragment->constrainedFlip,
          "fragment 5 flags");
    check(!atoms[8]->fragment->constrained &&
              !atoms[8]->fragment->constrainedFlip,
          "fragment 6 flags");
    for (auto* fragment : molecule->_fragments) {
        delete fragment;
    }
    molecule->_fragments.clear();
}

using Test = void (*)();

const std::pair<const char*, Test> tests[] = {
    {"Basics", Basics},
    {"Rings", Rings},
    {"Branching", Branching},
    {"BondOrder", BondOrder},
    {"testPolyominoCoordinatesOfSubstituent", testPolyominoCoordinatesOfSubstituent},
    {"testPolyominoSameAs", testPolyominoSameAs},
    {"testClockwiseOrderedSubstituents", testClockwiseOrderedSubstituents},
    {"testClockwiseOrderedNaN", testClockwiseOrderedNaN},
    {"testbicyclopentane", testbicyclopentane},
    {"testFusedRings", testFusedRings},
    {"testTemplates", testTemplates},
    {"testRingComplex", testRingComplex},
    {"testCoordgenFragmenter", testCoordgenFragmenter},
};

} // namespace

int main()
{
    for (const auto& test : tests) {
        try {
            test.second();
        } catch (const std::exception& error) {
            std::fprintf(stderr, "%s: %s\n", test.first, error.what());
            return 1;
        }
    }
    return 0;
}
