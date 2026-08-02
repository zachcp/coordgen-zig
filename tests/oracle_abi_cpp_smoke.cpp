#include "coordgen_abi.h"

#include <cstdint>

int main() {
    coordgen_atom_input_t atoms[2] = {};
    coordgen_bond_input_t bond = {};
    coordgen_input_t input = {};
    coordgen_result_t result = {};

    for (auto& atom : atoms) {
        atom.atomic_number = 6;
        atom.stereo_looking_from = COORDGEN_INVALID_INDEX;
        atom.stereo_atom_a = COORDGEN_INVALID_INDEX;
        atom.stereo_atom_b = COORDGEN_INVALID_INDEX;
    }
    bond.start = 0;
    bond.end = 1;
    bond.order = COORDGEN_BOND_SINGLE;
    bond.crossing_penalty_multiplier = 1.0f;
    input.options = coordgen_default_options();
    input.atoms = { atoms, 2, 0 };
    input.bonds = { &bond, 1, 0 };
    if (coordgen_generate(&input, &result) != COORDGEN_OK) return 1;
    const bool valid = result.coordinates.len == 2 &&
                       result.coordinates.ptr != nullptr &&
                       result.effective_bond_orders.len == 1 &&
                       result.effective_bond_orders.ptr[0] == COORDGEN_BOND_SINGLE;
    coordgen_result_free(&result);
    return valid && result.owner == nullptr ? 0 : 2;
}
