/*
 * Compiled and run by tools/check-install-isolation against a *scratch*
 * install prefix (`-I prefix/include -L prefix/lib -lcoordgen`), never
 * against this repository's build tree or source paths. Its only job is to
 * prove the installed artifacts are self-sufficient and behave as
 * documented for a consumer who has nothing but `$prefix`.
 */
#include "coordgen_abi.h"

#include <stddef.h>

int main(void) {
    coordgen_atom_input_t atoms[2] = {0};
    coordgen_bond_input_t bond = {0};
    coordgen_input_t input = {0};
    coordgen_result_t result = {0};

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
    input.atoms.ptr = atoms;
    input.atoms.len = 2;
    input.bonds.ptr = &bond;
    input.bonds.len = 1;

    if (coordgen_generate(&input, &result) != COORDGEN_OK) return 1;
    if (result.owner == NULL || result.coordinates.ptr == NULL ||
        result.coordinates.len != 2 || result.effective_bond_orders.len != 1) return 2;
    if (result.effective_bond_orders.ptr[0] != COORDGEN_BOND_SINGLE) return 3;
    coordgen_result_free(&result);
    if (result.owner != NULL || result.coordinates.ptr != NULL ||
        result.coordinates.len != 0 || result.effective_bond_orders.ptr != NULL ||
        result.effective_bond_orders.len != 0) return 4;

    return 0;
}
