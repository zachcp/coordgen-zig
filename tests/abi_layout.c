#include "coordgen_abi.h"

#include <stddef.h>

_Static_assert(sizeof(float) == 4, "CoordGen requires IEEE-width f32");
_Static_assert(sizeof(void *) == 8, "initial ABI targets are 64-bit");
_Static_assert(sizeof(coordgen_vec2_t) == 8, "coordgen_vec2_t layout");
_Static_assert(sizeof(coordgen_vec3_t) == 12, "coordgen_vec3_t layout");
_Static_assert(sizeof(coordgen_string_view_t) == 16, "string view layout");
_Static_assert(sizeof(coordgen_atom_input_t) == 52, "atom input layout");
_Static_assert(_Alignof(coordgen_atom_input_t) == 4, "atom input alignment");
_Static_assert(offsetof(coordgen_atom_input_t, template_coordinates) == 28,
               "atom template coordinate offset");
_Static_assert(offsetof(coordgen_atom_input_t, reserved) == 48,
               "atom reserved offset");
_Static_assert(sizeof(coordgen_bond_input_t) == 40, "bond input layout");
_Static_assert(sizeof(coordgen_residue_input_t) == 32, "residue input layout");
_Static_assert(sizeof(coordgen_residue_interaction_input_t) == 56,
               "residue interaction layout");
_Static_assert(sizeof(coordgen_options_t) == 56, "options layout");
_Static_assert(offsetof(coordgen_options_t, template_directory) == 40,
               "template directory offset");
_Static_assert(sizeof(coordgen_input_t) == 136, "input layout");
_Static_assert(_Alignof(coordgen_input_t) == 8, "input alignment");
_Static_assert(sizeof(coordgen_result_t) == 112, "result layout");
_Static_assert(offsetof(coordgen_result_t, owner) == 104, "result owner offset");
_Static_assert(COORDGEN_BOND_ZERO == 0, "zero-order bond value");
_Static_assert((int)COORDGEN_BOND_LENGTH == 50, "public bond length");

static int check_default_options(void) {
    coordgen_options_t options = coordgen_default_options();
    return options.precision == COORDGEN_PRECISION_STANDARD ? 0 : 1;
}

/*
 * The static asserts above check the header against itself. Calling into the
 * linked library is the other half: it turns a drift between this header's
 * declared coordgen_generate/coordgen_result_free signatures and what the
 * Zig implementation actually exports into a link failure, and it exercises
 * the documented error-code and ownership contracts end to end through the
 * real ABI boundary rather than by inspection.
 *
 * The discrete/continuous coordinate-generation pipeline is not wired yet
 * (see cgz-7v2.4), so the third case below - a structurally valid input -
 * deliberately expects COORDGEN_ERROR_UNSUPPORTED rather than coordinates.
 */
static int check_generate_error_paths(void) {
    coordgen_result_t result;

    coordgen_input_t empty_input = {0};
    empty_input.options = coordgen_default_options();
    if (coordgen_generate(&empty_input, &result) != COORDGEN_ERROR_EMPTY_GRAPH) return 1;
    /* "On failure, result is zeroed and requires no cleanup" (coordgen_abi.h) -
     * verify the zeroing, then verify that freeing it anyway is still safe. */
    if (result.owner != NULL || result.coordinates.ptr != NULL) return 1;
    coordgen_result_free(&result);

    coordgen_atom_input_t bad_atom = {0};
    bad_atom.atomic_number = 0; /* 0 is the internal virtual atom, rejected publicly */
    coordgen_input_t bad_atom_input = {0};
    bad_atom_input.options = coordgen_default_options();
    bad_atom_input.atoms.ptr = &bad_atom;
    bad_atom_input.atoms.len = 1;
    if (coordgen_generate(&bad_atom_input, &result) != COORDGEN_ERROR_INVALID_ATOMIC_NUMBER) {
        return 1;
    }
    coordgen_result_free(&result);

    coordgen_atom_input_t atoms[2] = {0};
    atoms[0].atomic_number = 6;
    atoms[0].stereo_looking_from = COORDGEN_INVALID_INDEX;
    atoms[0].stereo_atom_a = COORDGEN_INVALID_INDEX;
    atoms[0].stereo_atom_b = COORDGEN_INVALID_INDEX;
    atoms[1].atomic_number = 6;
    atoms[1].stereo_looking_from = COORDGEN_INVALID_INDEX;
    atoms[1].stereo_atom_a = COORDGEN_INVALID_INDEX;
    atoms[1].stereo_atom_b = COORDGEN_INVALID_INDEX;
    coordgen_bond_input_t bond = {0};
    bond.start = 0;
    bond.end = 1;
    bond.order = COORDGEN_BOND_SINGLE;
    bond.stereo_atom_a = COORDGEN_INVALID_INDEX;
    bond.stereo_atom_b = COORDGEN_INVALID_INDEX;
    coordgen_input_t valid_input = {0};
    valid_input.options = coordgen_default_options();
    valid_input.atoms.ptr = atoms;
    valid_input.atoms.len = 2;
    valid_input.bonds.ptr = &bond;
    valid_input.bonds.len = 1;
    if (coordgen_generate(&valid_input, &result) != COORDGEN_ERROR_UNSUPPORTED) return 1;
    coordgen_result_free(&result);

    return 0;
}

int main(void) {
    if (check_default_options() != 0) return 1;
    if (check_generate_error_paths() != 0) return 1;
    return 0;
}
