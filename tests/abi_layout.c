#include "coordgen_abi.h"

#include <stddef.h>

#define ASSERT_LAYOUT(type, expected_size, expected_align)                    \
    _Static_assert(sizeof(type) == (expected_size), #type " size");          \
    _Static_assert(_Alignof(type) == (expected_align), #type " alignment")
#define ASSERT_OFFSET(type, field, expected_offset)                           \
    _Static_assert(offsetof(type, field) == (expected_offset),                \
                   #type "." #field " offset")

_Static_assert(sizeof(float) == 4, "CoordGen requires IEEE-width f32");
_Static_assert(sizeof(void *) == 8, "initial ABI targets are 64-bit");

ASSERT_LAYOUT(coordgen_vec2_t, 8, 4);
ASSERT_OFFSET(coordgen_vec2_t, x, 0);
ASSERT_OFFSET(coordgen_vec2_t, y, 4);
ASSERT_LAYOUT(coordgen_vec3_t, 12, 4);
ASSERT_OFFSET(coordgen_vec3_t, x, 0);
ASSERT_OFFSET(coordgen_vec3_t, y, 4);
ASSERT_OFFSET(coordgen_vec3_t, z, 8);

ASSERT_LAYOUT(coordgen_string_view_t, 16, 8);
ASSERT_OFFSET(coordgen_string_view_t, ptr, 0);
ASSERT_OFFSET(coordgen_string_view_t, len, 8);
ASSERT_OFFSET(coordgen_string_view_t, reserved, 12);
ASSERT_LAYOUT(coordgen_index_span_t, 16, 8);
ASSERT_OFFSET(coordgen_index_span_t, ptr, 0);
ASSERT_OFFSET(coordgen_index_span_t, len, 8);
ASSERT_OFFSET(coordgen_index_span_t, reserved, 12);

ASSERT_LAYOUT(coordgen_atom_input_t, 52, 4);
ASSERT_OFFSET(coordgen_atom_input_t, atomic_number, 0);
ASSERT_OFFSET(coordgen_atom_input_t, formal_charge, 4);
ASSERT_OFFSET(coordgen_atom_input_t, flags, 8);
ASSERT_OFFSET(coordgen_atom_input_t, stereo, 12);
ASSERT_OFFSET(coordgen_atom_input_t, stereo_looking_from, 16);
ASSERT_OFFSET(coordgen_atom_input_t, stereo_atom_a, 20);
ASSERT_OFFSET(coordgen_atom_input_t, stereo_atom_b, 24);
ASSERT_OFFSET(coordgen_atom_input_t, template_coordinates, 28);
ASSERT_OFFSET(coordgen_atom_input_t, coordinates_3d, 36);
ASSERT_OFFSET(coordgen_atom_input_t, reserved, 48);

ASSERT_LAYOUT(coordgen_bond_input_t, 40, 4);
ASSERT_OFFSET(coordgen_bond_input_t, start, 0);
ASSERT_OFFSET(coordgen_bond_input_t, end, 4);
ASSERT_OFFSET(coordgen_bond_input_t, order, 8);
ASSERT_OFFSET(coordgen_bond_input_t, flags, 12);
ASSERT_OFFSET(coordgen_bond_input_t, stereo, 16);
ASSERT_OFFSET(coordgen_bond_input_t, stereo_atom_a, 20);
ASSERT_OFFSET(coordgen_bond_input_t, stereo_atom_b, 24);
ASSERT_OFFSET(coordgen_bond_input_t, display, 28);
ASSERT_OFFSET(coordgen_bond_input_t, crossing_penalty_multiplier, 32);
ASSERT_OFFSET(coordgen_bond_input_t, reserved, 36);

ASSERT_LAYOUT(coordgen_residue_input_t, 32, 8);
ASSERT_OFFSET(coordgen_residue_input_t, atom, 0);
ASSERT_OFFSET(coordgen_residue_input_t, residue_number, 4);
ASSERT_OFFSET(coordgen_residue_input_t, closest_ligand_atom, 8);
ASSERT_OFFSET(coordgen_residue_input_t, reserved, 12);
ASSERT_OFFSET(coordgen_residue_input_t, chain, 16);

ASSERT_LAYOUT(coordgen_residue_interaction_input_t, 56, 8);
ASSERT_OFFSET(coordgen_residue_interaction_input_t, start, 0);
ASSERT_OFFSET(coordgen_residue_interaction_input_t, end, 4);
ASSERT_OFFSET(coordgen_residue_interaction_input_t, other_start_atoms, 8);
ASSERT_OFFSET(coordgen_residue_interaction_input_t, other_end_atoms, 24);
ASSERT_OFFSET(coordgen_residue_interaction_input_t, crossing_penalty_multiplier, 40);
ASSERT_OFFSET(coordgen_residue_interaction_input_t, reserved, 44);

#define ASSERT_SPAN(type)                                                     \
    ASSERT_LAYOUT(type, 16, 8);                                               \
    ASSERT_OFFSET(type, ptr, 0);                                              \
    ASSERT_OFFSET(type, len, 8);                                              \
    ASSERT_OFFSET(type, reserved, 12)
ASSERT_SPAN(coordgen_atom_span_t);
ASSERT_SPAN(coordgen_bond_span_t);
ASSERT_SPAN(coordgen_residue_span_t);
ASSERT_SPAN(coordgen_residue_interaction_span_t);
ASSERT_SPAN(coordgen_vec2_span_t);
ASSERT_SPAN(coordgen_u32_span_t);
#undef ASSERT_SPAN

ASSERT_LAYOUT(coordgen_options_t, 56, 8);
ASSERT_OFFSET(coordgen_options_t, precision, 0);
ASSERT_OFFSET(coordgen_options_t, score_residue_interactions, 4);
ASSERT_OFFSET(coordgen_options_t, treat_nonterminal_bonds_to_metal_as_zero_order, 8);
ASSERT_OFFSET(coordgen_options_t, even_angles, 12);
ASSERT_OFFSET(coordgen_options_t, skip_minimization, 16);
ASSERT_OFFSET(coordgen_options_t, force_open_macrocycles, 20);
ASSERT_OFFSET(coordgen_options_t, constrain_all_atoms, 24);
/* cgz-r25 renamed this slot from build_from_fragments to reserved. The
 * offset is unchanged, which is the whole point of asserting it here: the
 * rename is a name-only amendment to the frozen layout. */
ASSERT_OFFSET(coordgen_options_t, reserved, 28);
ASSERT_OFFSET(coordgen_options_t, debug_coordinates, 32);
ASSERT_OFFSET(coordgen_options_t, load_templates, 36);
ASSERT_OFFSET(coordgen_options_t, template_directory, 40);

ASSERT_LAYOUT(coordgen_input_t, 136, 8);
ASSERT_OFFSET(coordgen_input_t, options, 0);
ASSERT_OFFSET(coordgen_input_t, atoms, 56);
ASSERT_OFFSET(coordgen_input_t, bonds, 72);
ASSERT_OFFSET(coordgen_input_t, residues, 88);
ASSERT_OFFSET(coordgen_input_t, residue_interactions, 104);
ASSERT_OFFSET(coordgen_input_t, extra_bonds, 120);

ASSERT_LAYOUT(coordgen_result_t, 112, 8);
ASSERT_OFFSET(coordgen_result_t, coordinates, 0);
ASSERT_OFFSET(coordgen_result_t, input_to_internal, 16);
ASSERT_OFFSET(coordgen_result_t, internal_to_input, 32);
ASSERT_OFFSET(coordgen_result_t, effective_bond_orders, 48);
ASSERT_OFFSET(coordgen_result_t, bond_displays, 64);
ASSERT_OFFSET(coordgen_result_t, atom_stereo, 80);
ASSERT_OFFSET(coordgen_result_t, clean_pose, 96);
ASSERT_OFFSET(coordgen_result_t, reserved, 100);
ASSERT_OFFSET(coordgen_result_t, owner, 104);

#undef ASSERT_OFFSET
#undef ASSERT_LAYOUT
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
