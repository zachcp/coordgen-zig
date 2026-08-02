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

int main(void) {
    coordgen_options_t options = coordgen_default_options();
    return options.precision == COORDGEN_PRECISION_STANDARD ? 0 : 1;
}
