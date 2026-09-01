/*
 * Out-of-tree C++ consumer for the installed C ABI. This deliberately uses
 * only the installed C header and archive; the optional C++ facade does not
 * exist yet, so C++ interoperability means the stable extern "C" surface is
 * usable with ordinary C++17 types and initialization.
 */
#include "coordgen_abi.h"

#include <type_traits>

static_assert(std::is_standard_layout<coordgen_input_t>::value,
              "input DTO must remain standard-layout in C++");
static_assert(std::is_standard_layout<coordgen_result_t>::value,
              "result DTO must remain standard-layout in C++");

int main() {
    coordgen_atom_input_t atoms[2] = {};
    coordgen_bond_input_t bond = {};
    coordgen_input_t input = {};
    coordgen_result_t result = {};

    atoms[0].atomic_number = 6;
    atoms[0].stereo_looking_from = COORDGEN_INVALID_INDEX;
    atoms[0].stereo_atom_a = COORDGEN_INVALID_INDEX;
    atoms[0].stereo_atom_b = COORDGEN_INVALID_INDEX;
    atoms[1] = atoms[0];
    bond.start = 0;
    bond.end = 1;
    bond.order = COORDGEN_BOND_SINGLE;
    bond.stereo_atom_a = COORDGEN_INVALID_INDEX;
    bond.stereo_atom_b = COORDGEN_INVALID_INDEX;

    input.options = coordgen_default_options();
    input.atoms = {atoms, 2, 0};
    input.bonds = {&bond, 1, 0};

    if (input.options.precision != COORDGEN_PRECISION_STANDARD) return 1;
    if (coordgen_generate(&input, &result) != COORDGEN_OK) return 2;
    if (result.owner == nullptr || result.coordinates.ptr == nullptr ||
        result.coordinates.len != 2 || result.effective_bond_orders.len != 1) return 3;
    if (result.effective_bond_orders.ptr[0] != COORDGEN_BOND_SINGLE) return 4;

    coordgen_result_free(&result);
    return result.owner == nullptr && result.coordinates.ptr == nullptr &&
                   result.coordinates.len == 0 &&
                   result.effective_bond_orders.ptr == nullptr &&
                   result.effective_bond_orders.len == 0
               ? 0
               : 5;
}
