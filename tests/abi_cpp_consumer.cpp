#include "coordgen_abi.h"

#include <type_traits>

static_assert(std::is_standard_layout<coordgen_input_t>::value,
              "C++ consumer requires standard-layout DTOs");
static_assert(std::is_standard_layout<coordgen_result_t>::value,
              "C++ consumer requires standard-layout results");

/*
 * Compile-time shape of the dominant RDKit flow: atom flags/template data,
 * bonds and stereo are present before the one-shot call; output is addressed
 * by the original atom array position and scaled by the public bond length.
 * Linking/running this against the pinned C++ oracle belongs to cgz-r04/r03.
 */
void rdkit_shaped_consumer(void) {
    coordgen_atom_input_t atoms[4] = {};
    for (coordgen_atom_input_t &atom : atoms) {
        atom.atomic_number = 6;
        atom.stereo_looking_from = COORDGEN_INVALID_INDEX;
        atom.stereo_atom_a = COORDGEN_INVALID_INDEX;
        atom.stereo_atom_b = COORDGEN_INVALID_INDEX;
    }
    atoms[0].formal_charge = -1;
    atoms[0].flags = COORDGEN_ATOM_FIXED |
                     COORDGEN_ATOM_CONSTRAINED |
                     COORDGEN_ATOM_HAS_TEMPLATE_COORDINATES;
    atoms[0].template_coordinates = {10.0f, 20.0f};

    coordgen_bond_input_t bond = {};
    bond.start = 0;
    bond.end = 1;
    bond.order = COORDGEN_BOND_DOUBLE;
    bond.stereo = COORDGEN_BOND_STEREO_CIS;
    bond.stereo_atom_a = 2;
    bond.stereo_atom_b = 3;
    bond.crossing_penalty_multiplier = 1.0f;

    coordgen_input_t input = {};
    input.options = coordgen_default_options();
    input.options.treat_nonterminal_bonds_to_metal_as_zero_order = 1;
    input.atoms = {atoms, 4, 0};
    input.bonds = {&bond, 1, 0};

    coordgen_result_t result = {};
    coordgen_error_t (*generate)(const coordgen_input_t *, coordgen_result_t *) =
        &coordgen_generate;
    void (*free_result)(coordgen_result_t *) = &coordgen_result_free;
    (void)generate;
    (void)free_result;
    (void)result;
    (void)COORDGEN_BOND_LENGTH;
}
