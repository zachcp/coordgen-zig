/*
 * Fixture-backed half of the Boost/maeparser-free public-test rehost. The
 * arguments are normalized fixture dumps produced by the repo-owned Zig MAE
 * reader, in this fixed order: sample, templates, chirality, nonterminal
 * metal, terminal metal, macrocycle, and test_mol.
 */

#include <cmath>
#include <cstdio>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

#include "CoordgenFragmentBuilder.h"
#include "CoordgenMacrocycleBuilder.h"
#include "sketcherMinimizer.h"
#include "sketcherMinimizerAtom.h"
#include "sketcherMinimizerBendInteraction.h"
#include "sketcherMinimizerBond.h"
#include "sketcherMinimizerMolecule.h"
#include "sketcherMinimizerMaths.h"
#include "sketcherMinimizerStretchInteraction.h"

namespace
{

struct Atom {
    unsigned atomic_number;
    float x;
    float y;
};

struct Bond {
    unsigned from;
    unsigned to;
    unsigned order;
};

struct Structure {
    std::vector<Atom> atoms;
    std::vector<Bond> bonds;
};

void check(bool condition, const char* message)
{
    if (!condition) {
        throw std::runtime_error(message);
    }
}

float fromBits(std::uint32_t bits)
{
    float value;
    static_assert(sizeof(value) == sizeof(bits), "f32 width");
    std::memcpy(&value, &bits, sizeof(value));
    return value;
}

std::vector<Structure> readStructures(const char* path)
{
    std::ifstream input(path);
    check(input.good(), "cannot open normalized fixture");
    std::vector<Structure> structures;
    std::string tag;
    while (input >> tag) {
        check(tag == "structure", "expected structure record");
        std::size_t atom_count = 0;
        std::size_t bond_count = 0;
        input >> atom_count >> bond_count;
        Structure structure;
        structure.atoms.reserve(atom_count);
        structure.bonds.reserve(bond_count);
        for (std::size_t index = 0; index < atom_count; ++index) {
            unsigned atomic_number = 0;
            std::uint32_t x = 0;
            std::uint32_t y = 0;
            input >> tag >> atomic_number >> x >> y;
            check(tag == "atom", "expected atom record");
            structure.atoms.push_back({atomic_number, fromBits(x), fromBits(y)});
        }
        for (std::size_t index = 0; index < bond_count; ++index) {
            Bond bond{};
            input >> tag >> bond.from >> bond.to >> bond.order;
            check(tag == "bond", "expected bond record");
            check(bond.from < atom_count && bond.to < atom_count,
                  "normalized bond index out of range");
            structure.bonds.push_back(bond);
        }
        input >> tag;
        check(tag == "end", "expected end record");
        structures.push_back(std::move(structure));
    }
    return structures;
}

sketcherMinimizerMolecule* makeMolecule(const Structure& structure)
{
    auto* molecule = new sketcherMinimizerMolecule();
    for (const auto& source : structure.atoms) {
        auto* atom = molecule->addNewAtom();
        atom->setAtomicNumber(static_cast<int>(source.atomic_number));
        atom->setCoordinates({source.x, source.y});
    }
    for (const auto& source : structure.bonds) {
        auto* bond = molecule->addNewBond(
            molecule->getAtoms().at(source.from),
            molecule->getAtoms().at(source.to));
        bond->setBondOrder(static_cast<int>(source.order));
    }
    return molecule;
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
                    std::map<sketcherMinimizerAtom*, int>& indices)
{
    const float target = BONDLENGTH * BONDLENGTH;
    const float tolerance = target * 0.1f;
    for (auto* bond : molecule.getBonds()) {
        const float distance = sketcherMinimizerMaths::squaredDistance(
            bond->getStartAtom()->getCoordinates(),
            bond->getEndAtom()->getCoordinates());
        if (distance < target - tolerance || distance > target + tolerance) {
            (void)indices;
            return false;
        }
    }
    return true;
}

void SampleTest(const std::vector<Structure>& structures)
{
    check(structures.size() == 1, "SampleTest structure count");
    auto* molecule = makeMolecule(structures[0]);
    check(molecule->getAtoms().size() == 26, "SampleTest atom count");
    check(molecule->getBonds().size() == 26, "SampleTest bond count");
    sketcherMinimizer minimizer;
    minimizer.initialize(molecule);
    minimizer.runGenerateCoordinates();
    for (auto* atom : molecule->getAtoms()) {
        const auto coordinates = atom->getCoordinates();
        check(coordinates.x() != 0 || coordinates.y() != 0,
              "SampleTest generated coordinate");
    }
    auto indices = reportingIndices(*molecule);
    check(bondsNearIdeal(*molecule, indices), "SampleTest bond lengths");
}

void TemplateTest(const std::vector<Structure>& structures)
{
    check(structures.size() == 82, "TemplateTest structure count");
    const std::unordered_set<std::size_t> no_match{
        1, 8, 19, 20, 22, 32, 43, 53, 65, 66, 67,
    };
    const std::unordered_set<std::size_t> match_incorrectly{18, 27};
    for (std::size_t index = 0; index < structures.size(); ++index) {
        if (no_match.count(index) != 0) {
            continue;
        }
        auto* molecule = makeMolecule(structures[index]);
        const std::size_t atom_count = molecule->getAtoms().size();
        sketcherMinimizer minimizer;
        minimizer.initialize(molecule);
        minimizer.runGenerateCoordinates();
        check(atom_count == molecule->getAtoms().size(),
              "TemplateTest preserves atom count");
        bool any_rigid = false;
        bool all_rigid = true;
        for (auto* atom : molecule->getAtoms()) {
            any_rigid = any_rigid || atom->rigid;
            all_rigid = all_rigid && atom->rigid;
        }
        check(any_rigid, "TemplateTest finds template");
        if (match_incorrectly.count(index) == 0) {
            check(all_rigid, "TemplateTest makes all atoms rigid");
        }
    }
}

void ClearWedgesTest(const std::vector<Structure>& structures)
{
    check(structures.size() == 1, "ClearWedgesTest structure count");
    auto* molecule = makeMolecule(structures[0]);
    check(molecule->getBonds().size() >= 4, "ClearWedgesTest bond count");
    for (std::size_t index = 0; index < 4; ++index) {
        molecule->getBonds().at(index)->hasStereochemistryDisplay = true;
    }
    auto* carbon = molecule->getAtoms().at(0);
    check(carbon->atomicNumber == 6, "ClearWedgesTest carbon");
    carbon->hasStereochemistrySet = true;
    carbon->isR = true;
    sketcherMinimizer minimizer;
    minimizer.initialize(molecule);
    minimizer.runGenerateCoordinates();
    check(!molecule->getBonds().at(0)->hasStereochemistryDisplay,
          "ClearWedgesTest clears first wedge");
    check(molecule->getBonds().at(1)->hasStereochemistryDisplay,
          "ClearWedgesTest writes second wedge");
    check(!molecule->getBonds().at(2)->hasStereochemistryDisplay,
          "ClearWedgesTest clears third wedge");
    check(molecule->getBonds().at(3)->hasStereochemistryDisplay,
          "ClearWedgesTest writes fourth wedge");
}

void checkMetal(const std::vector<Structure>& structures, bool automatic_zero_order)
{
    check(structures.size() == 1, "metal structure count");
    auto* molecule = makeMolecule(structures[0]);
    auto* aluminum = molecule->getAtoms().at(0);
    auto* nitrogen = molecule->getAtoms().at(1);
    check(aluminum->atomicNumber == 13, "metal aluminum");
    check(nitrogen->atomicNumber == 7, "metal nitrogen");
    sketcherMinimizer minimizer;
    if (!automatic_zero_order) {
        minimizer.setTreatNonterminalBondsToMetalAsZOBs(false);
    }
    minimizer.initialize(molecule);
    minimizer.runGenerateCoordinates();
    const float length = (aluminum->coordinates - nitrogen->coordinates).length();
    check(length > 48 && length < 52, "metal bond length");
    auto indices = reportingIndices(*molecule);
    check(bondsNearIdeal(*molecule, indices), "metal ideal bonds");
}

void DisableMetalZOBs(const std::vector<Structure>& structures)
{
    checkMetal(structures, false);
}

void terminalMetalZOBs(const std::vector<Structure>& structures)
{
    checkMetal(structures, true);
}

void testMinimizedRingsShape(const std::vector<Structure>& structures)
{
    check(structures.size() == 1, "ring shape structure count");
    auto* molecule = makeMolecule(structures[0]);
    sketcherMinimizer minimizer;
    minimizer.initialize(molecule);
    minimizer.runGenerateCoordinates();
    for (auto* interaction : minimizer.getStretchInteractions()) {
        auto* ring = sketcherMinimizer::sameRing(interaction->atom1,
                                                  interaction->atom2);
        if (ring != nullptr && !ring->isMacrocycle()) {
            const float length = (interaction->atom1->coordinates -
                                  interaction->atom2->coordinates).length();
            check(length > 48 && length < 52, "ring stretch length");
        }
    }
    int angle_count = 0;
    for (auto* interaction : minimizer.getBendInteractions()) {
        if (!interaction->isRing) {
            continue;
        }
        auto* ring = sketcherMinimizer::sameRing(
            interaction->atom1, interaction->atom2, interaction->atom3);
        check(ring != nullptr, "ring bend membership");
        check(!ring->isMacrocycle(), "ring bend non-macrocycle");
        const float angle = interaction->angle();
        check(angle > interaction->restV - 2 && angle < interaction->restV + 2,
              "ring bend angle");
        ++angle_count;
    }
    check(angle_count == 32, "ring bend count");
}

void testGetDoubleBondConstraints(const std::vector<Structure>& structures)
{
    check(structures.size() == 1, "double-bond structure count");
    auto* molecule = makeMolecule(structures[0]);
    sketcherMinimizer minimizer;
    CoordgenFragmentBuilder fragment_builder;
    CoordgenMacrocycleBuilder macrocycle_builder;
    minimizer.initialize(molecule);
    for (auto* component : minimizer.getMolecules()) {
        for (auto* ring : component->getRings()) {
            auto atoms = fragment_builder.orderRingAtoms(ring);
            const auto constraints =
                macrocycle_builder.getDoubleBondConstraints(atoms);
            check(constraints.empty(), "small-ring double-bond constraints");
        }
    }
}

} // namespace

int main(int argc, char** argv)
{
    if (argc != 8) {
        return 2;
    }
    try {
        SampleTest(readStructures(argv[1]));
        TemplateTest(readStructures(argv[2]));
        ClearWedgesTest(readStructures(argv[3]));
        DisableMetalZOBs(readStructures(argv[4]));
        terminalMetalZOBs(readStructures(argv[5]));
        testMinimizedRingsShape(readStructures(argv[6]));
        testGetDoubleBondConstraints(readStructures(argv[7]));
    } catch (const std::exception& error) {
        std::fprintf(stderr, "%s\n", error.what());
        return 1;
    }
    return 0;
}
