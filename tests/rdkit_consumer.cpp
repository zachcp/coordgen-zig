/*
 * cgz-r07 sufficiency test: the dominant downstream consumer's flow, executed
 * against the pinned C++ facade through the stable C ABI.
 *
 * RDKit's External/CoordGen/CoordGen.h builds atoms and bonds, sets per-atom
 * flags and template coordinates, assigns bond stereo, runs the generator, and
 * then reads coordinates back *from its own atom pointers* - that is, in its
 * own input order - dividing by a bond length of 50. Pointers made that order
 * contract invisible upstream. An index-based ABI has to state it, and this
 * test is what makes the statement falsifiable: the atoms are handed over in an
 * order that is deliberately not the connectivity order, so if the result were
 * ever returned in canonical or component order the bonded-distance checks
 * below would fail rather than silently produce transposed coordinates.
 *
 * tests/abi_cpp_consumer.cpp covers the same surface at compile time only, for
 * builds that must not depend on the oracle.
 */

#include "coordgen_abi.h"

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <vector>

namespace {

int failures = 0;

void check(bool condition, const char *what) {
    if (!condition) {
        std::fprintf(stderr, "rdkit-consumer: FAILED %s\n", what);
        ++failures;
    }
}

void checkClose(float actual, float expected, float tolerance, const char *what) {
    if (!(std::fabs(actual - expected) <= tolerance)) {
        std::fprintf(stderr, "rdkit-consumer: FAILED %s (got %.6f, expected %.6f +/- %.6f)\n",
                     what, static_cast<double>(actual), static_cast<double>(expected),
                     static_cast<double>(tolerance));
        ++failures;
    }
}

float distance(const coordgen_vec2_t &a, const coordgen_vec2_t &b) {
    const float dx = a.x - b.x;
    const float dy = a.y - b.y;
    return std::sqrt(dx * dx + dy * dy);
}

coordgen_atom_input_t carbon() {
    coordgen_atom_input_t atom = {};
    atom.atomic_number = 6;
    atom.stereo_looking_from = COORDGEN_INVALID_INDEX;
    atom.stereo_atom_a = COORDGEN_INVALID_INDEX;
    atom.stereo_atom_b = COORDGEN_INVALID_INDEX;
    return atom;
}

coordgen_bond_input_t bond(uint32_t start, uint32_t end, coordgen_bond_order_t order) {
    coordgen_bond_input_t value = {};
    value.start = start;
    value.end = end;
    value.order = order;
    value.stereo_atom_a = COORDGEN_INVALID_INDEX;
    value.stereo_atom_b = COORDGEN_INVALID_INDEX;
    value.crossing_penalty_multiplier = 1.0f;
    return value;
}

/*
 * Toluene, kekulized, with the ring closed through atoms the caller lists out
 * of order: input index 0 is the methyl carbon, and the ring is enumerated as
 * 3, 6, 1, 5, 2, 4. No consecutive pair of input indices is a bond except by
 * accident, so input order and any internal ordering are distinguishable.
 */
struct Molecule {
    std::vector<coordgen_atom_input_t> atoms;
    std::vector<coordgen_bond_input_t> bonds;
};

Molecule toluene() {
    Molecule molecule;
    molecule.atoms.assign(7, carbon());

    // Ring, in caller order: 3-6, 6-1, 1-5, 5-2, 2-4, 4-3. Methyl: 0-3.
    molecule.bonds.push_back(bond(3, 6, COORDGEN_BOND_DOUBLE));
    molecule.bonds.push_back(bond(6, 1, COORDGEN_BOND_SINGLE));
    molecule.bonds.push_back(bond(1, 5, COORDGEN_BOND_DOUBLE));
    molecule.bonds.push_back(bond(5, 2, COORDGEN_BOND_SINGLE));
    molecule.bonds.push_back(bond(2, 4, COORDGEN_BOND_DOUBLE));
    molecule.bonds.push_back(bond(4, 3, COORDGEN_BOND_SINGLE));
    molecule.bonds.push_back(bond(0, 3, COORDGEN_BOND_SINGLE));
    return molecule;
}

coordgen_input_t inputFor(const Molecule &molecule) {
    coordgen_input_t input = {};
    input.options = coordgen_default_options();
    // The two options RDKit actually sets.
    input.options.treat_nonterminal_bonds_to_metal_as_zero_order = 1;
    input.options.precision = COORDGEN_PRECISION_STANDARD;
    input.atoms = { molecule.atoms.data(), static_cast<uint32_t>(molecule.atoms.size()), 0 };
    input.bonds = { molecule.bonds.data(), static_cast<uint32_t>(molecule.bonds.size()), 0 };
    return input;
}

/* Output is addressed by input atom index, and the permutation maps that make
 * the internal reordering observable are mutual inverses. */
void checkInputOrderContract(const Molecule &molecule, const coordgen_result_t &result) {
    const uint32_t atom_count = static_cast<uint32_t>(molecule.atoms.size());
    check(result.coordinates.len == atom_count, "one coordinate per input atom");
    check(result.input_to_internal.len == atom_count, "input_to_internal covers every input atom");
    check(result.internal_to_input.len == atom_count, "internal_to_input covers every atom");
    check(result.effective_bond_orders.len == molecule.bonds.size(),
          "one effective bond order per input bond");
    if (failures != 0) return;

    for (uint32_t input_index = 0; input_index < atom_count; ++input_index) {
        const uint32_t internal = result.input_to_internal.ptr[input_index];
        check(internal < atom_count, "input_to_internal stays in range");
        if (internal >= atom_count) return;
        check(result.internal_to_input.ptr[internal] == input_index,
              "internal_to_input inverts input_to_internal");
    }

    // Bonded pairs are one bond length apart *when addressed by input index*.
    // Returning any other order would break this for a scrambled input.
    // Measured against the pinned facade: 49.996 to 50.001 across these seven
    // bonds. The tolerance is an order of magnitude above that spread and three
    // below the value being asserted, so a transposed result cannot pass.
    for (const coordgen_bond_input_t &b : molecule.bonds) {
        const float length = distance(result.coordinates.ptr[b.start], result.coordinates.ptr[b.end]);
        checkClose(length, COORDGEN_BOND_LENGTH, 0.05f, "bonded atoms are one bond length apart");
    }

    // RDKit divides by 50 to reach unit bond length; the constant is public so
    // consumers do not have to rediscover it.
    const coordgen_bond_input_t &methyl = molecule.bonds.back();
    const float scaled = distance(result.coordinates.ptr[methyl.start],
                                  result.coordinates.ptr[methyl.end]) /
                         COORDGEN_BOND_LENGTH;
    checkClose(scaled, 1.0f, 0.001f, "scaling by COORDGEN_BOND_LENGTH yields unit bonds");

    // The methyl carbon is a substituent, not a ring member: it is farther from
    // the far side of the ring than any bonded pair. A transposed result would
    // not place it consistently.
    check(distance(result.coordinates.ptr[0], result.coordinates.ptr[5]) > COORDGEN_BOND_LENGTH,
          "non-bonded atoms are farther apart than bonded atoms");

    for (uint32_t i = 0; i < molecule.bonds.size(); ++i) {
        check(result.effective_bond_orders.ptr[i] == molecule.bonds[i].order,
              "effective bond orders keep input bond order and value for an all-carbon molecule");
    }
}

/* The whole input is one borrowed struct, so RDKit's build-bonds ->
 * assignBondsAndNeighbors -> set-stereo sequence has no representation here:
 * there is no order in which a caller can get it wrong. */
int generateOnce(const Molecule &molecule) {
    const coordgen_input_t input = inputFor(molecule);
    coordgen_result_t result = {};
    const coordgen_error_t status = coordgen_generate(&input, &result);
    if (status != COORDGEN_OK) {
        std::fprintf(stderr, "rdkit-consumer: coordgen_generate failed with %u\n", status);
        return 1;
    }
    checkInputOrderContract(molecule, result);
    coordgen_result_free(&result);
    check(result.owner == nullptr, "freeing a result clears its owner");
    check(result.coordinates.ptr == nullptr, "freeing a result clears its spans");
    return 0;
}

/* Same molecule, atoms renumbered by an unrelated permutation. Coordinates are
 * a rigid-motion equivalence class, so the invariant that must survive is the
 * geometry: every bond length and the full distance matrix, read through each
 * input's own indices. */
int checkOrderIndependence() {
    const Molecule original = toluene();
    const uint32_t permutation[7] = { 4, 0, 6, 2, 5, 1, 3 };

    Molecule permuted;
    permuted.atoms.assign(original.atoms.size(), carbon());
    for (const coordgen_bond_input_t &b : original.bonds) {
        permuted.bonds.push_back(bond(permutation[b.start], permutation[b.end], b.order));
    }

    const coordgen_input_t first_input = inputFor(original);
    const coordgen_input_t second_input = inputFor(permuted);
    coordgen_result_t first = {};
    coordgen_result_t second = {};
    if (coordgen_generate(&first_input, &first) != COORDGEN_OK) return 1;
    if (coordgen_generate(&second_input, &second) != COORDGEN_OK) {
        coordgen_result_free(&first);
        return 1;
    }

    for (uint32_t i = 0; i < original.atoms.size(); ++i) {
        for (uint32_t j = i + 1; j < original.atoms.size(); ++j) {
            const float lhs = distance(first.coordinates.ptr[i], first.coordinates.ptr[j]);
            const float rhs = distance(second.coordinates.ptr[permutation[i]],
                                       second.coordinates.ptr[permutation[j]]);
            // Measured deviation on this molecule is exactly 0. The tolerance
            // is kept nonzero because the general claim is geometric identity,
            // not bit identity: which side of the discrete search a molecule
            // lands on is what cgz-r13's parity ceiling is about.
            checkClose(rhs, lhs, 0.001f, "distances survive renumbering of the input");
        }
    }

    coordgen_result_free(&first);
    coordgen_result_free(&second);
    return 0;
}

} // namespace

int main() {
    const Molecule molecule = toluene();
    if (generateOnce(molecule) != 0) return 1;
    if (checkOrderIndependence() != 0) return 1;
    if (failures != 0) {
        std::fprintf(stderr, "rdkit-consumer: %d contract check(s) failed\n", failures);
        return 2;
    }
    return 0;
}
