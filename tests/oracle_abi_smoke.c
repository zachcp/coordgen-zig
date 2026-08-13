#include "coordgen_abi.h"
#include "coordgen_probe.h"

#include <stdint.h>

static int stable_api_smoke(void) {
    coordgen_atom_input_t atoms[3] = {0};
    coordgen_bond_input_t bonds[3] = {0};
    coordgen_input_t input = {0};
    coordgen_result_t result = {0};
    uint32_t index;

    for (index = 0; index < 3; ++index) {
        atoms[index].atomic_number = 6;
        atoms[index].stereo_looking_from = COORDGEN_INVALID_INDEX;
        atoms[index].stereo_atom_a = COORDGEN_INVALID_INDEX;
        atoms[index].stereo_atom_b = COORDGEN_INVALID_INDEX;
    }
    bonds[0] = (coordgen_bond_input_t){0, 1, COORDGEN_BOND_SINGLE, 0, 0,
                                        COORDGEN_INVALID_INDEX, COORDGEN_INVALID_INDEX, 0, 1.0f, 0};
    bonds[1] = (coordgen_bond_input_t){1, 2, COORDGEN_BOND_SINGLE, 0, 0,
                                        COORDGEN_INVALID_INDEX, COORDGEN_INVALID_INDEX, 0, 1.0f, 0};
    bonds[2] = (coordgen_bond_input_t){2, 0, COORDGEN_BOND_SINGLE, 0, 0,
                                        COORDGEN_INVALID_INDEX, COORDGEN_INVALID_INDEX, 0, 1.0f, 0};
    input.options = coordgen_default_options();
    input.atoms = (coordgen_atom_span_t){atoms, 3, 0};
    input.bonds = (coordgen_bond_span_t){bonds, 3, 0};
    if (coordgen_generate(&input, &result) != COORDGEN_OK) return 1;
    if (result.coordinates.len != 3 || result.input_to_internal.len != 3 ||
        result.internal_to_input.len != 3 || result.effective_bond_orders.len != 3 ||
        result.owner == 0) {
        coordgen_result_free(&result);
        return 2;
    }
    for (index = 0; index < result.internal_to_input.len; ++index) {
        if (result.internal_to_input.ptr[index] >= 3 ||
            result.input_to_internal.ptr[result.internal_to_input.ptr[index]] != index) {
            coordgen_result_free(&result);
            return 3;
        }
    }
    coordgen_result_free(&result);
    return result.owner == 0 ? 0 : 4;
}

static int probe_api_smoke(void) {
    coordgen_atom_input_t atoms[7] = {0};
    coordgen_bond_input_t bonds[7] = {0};
    coordgen_input_t input = {0};
    coordgen_probe_result_t probe = {0};
    uint32_t index;

    int found_template = 0;
    int moved_component = 0;
    for (index = 0; index < 7; ++index) {
        atoms[index].atomic_number = 6;
        atoms[index].stereo_looking_from = COORDGEN_INVALID_INDEX;
        atoms[index].stereo_atom_a = COORDGEN_INVALID_INDEX;
        atoms[index].stereo_atom_b = COORDGEN_INVALID_INDEX;
    }
    bonds[0] = (coordgen_bond_input_t){0, 1, COORDGEN_BOND_SINGLE, 0, 0,
                                        COORDGEN_INVALID_INDEX, COORDGEN_INVALID_INDEX, 0, 1.0f, 0};
    bonds[1] = (coordgen_bond_input_t){0, 3, COORDGEN_BOND_SINGLE, 0, 0,
                                        COORDGEN_INVALID_INDEX, COORDGEN_INVALID_INDEX, 0, 1.0f, 0};
    bonds[2] = (coordgen_bond_input_t){1, 2, COORDGEN_BOND_SINGLE, 0, 0,
                                        COORDGEN_INVALID_INDEX, COORDGEN_INVALID_INDEX, 0, 1.0f, 0};
    bonds[3] = (coordgen_bond_input_t){1, 4, COORDGEN_BOND_SINGLE, 0, 0,
                                        COORDGEN_INVALID_INDEX, COORDGEN_INVALID_INDEX, 0, 1.0f, 0};
    bonds[4] = (coordgen_bond_input_t){2, 3, COORDGEN_BOND_SINGLE, 0, 0,
                                        COORDGEN_INVALID_INDEX, COORDGEN_INVALID_INDEX, 0, 1.0f, 0};
    bonds[5] = (coordgen_bond_input_t){3, 4, COORDGEN_BOND_SINGLE, 0, 0,
                                        COORDGEN_INVALID_INDEX, COORDGEN_INVALID_INDEX, 0, 1.0f, 0};
    bonds[6] = (coordgen_bond_input_t){5, 6, COORDGEN_BOND_SINGLE, 0, 0,
                                        COORDGEN_INVALID_INDEX, COORDGEN_INVALID_INDEX, 0, 1.0f, 0};
    input.options = coordgen_default_options();
    input.atoms = (coordgen_atom_span_t){atoms, 7, 0};
    input.bonds = (coordgen_bond_span_t){bonds, 7, 0};
    if (coordgen_probe_generate(&input, &probe) != COORDGEN_OK) return 10;
    if (probe.input_to_internal.len != 7 || probe.internal_to_input.len != 7 ||
        probe.morgan_ranks.len != 7 || probe.ring_count == 0 ||
        probe.fragment_count == 0 || probe.component_count != 2 ||
        probe.component_atoms.len != 7 ||
        probe.owner == 0) {
        coordgen_probe_result_free(&probe);
        return 11;
    }
    for (index = 0; index < probe.fragment_count; ++index) {
        if (probe.fragments[index].template_match == 80 &&
            probe.fragments[index].template_mapping_count == 5 &&
            probe.template_mapping_count >= probe.fragments[index].template_mapping_count &&
            probe.template_mapping[probe.fragments[index].template_mapping_start].input_atom < 7) {
            found_template = 1;
        }
    }
    for (index = 0; index < probe.component_count; ++index) {
        if (probe.components[index].transform_status != COORDGEN_PROBE_TRANSFORM_OBSERVED) {
            coordgen_probe_result_free(&probe);
            return 12;
        }
        if (probe.components[index].transform[4] != 0.0f ||
            probe.components[index].transform[5] != 0.0f) {
            moved_component = 1;
        }
    }
    if (!found_template || !moved_component) {
        coordgen_probe_result_free(&probe);
        return 13;
    }
    coordgen_probe_result_free(&probe);
    return probe.owner == 0 ? 0 : 14;
}

static int no_template_probe_smoke(void) {
    coordgen_atom_input_t atoms[2] = {0};
    coordgen_bond_input_t bond = {0};
    coordgen_input_t input = {0};
    coordgen_probe_result_t probe = {0};
    uint32_t index;

    for (index = 0; index < 2; ++index) {
        atoms[index].atomic_number = 6;
        atoms[index].stereo_looking_from = COORDGEN_INVALID_INDEX;
        atoms[index].stereo_atom_a = COORDGEN_INVALID_INDEX;
        atoms[index].stereo_atom_b = COORDGEN_INVALID_INDEX;
    }
    bond.start = 0;
    bond.end = 1;
    bond.order = COORDGEN_BOND_SINGLE;
    bond.crossing_penalty_multiplier = 1.0f;
    input.options = coordgen_default_options();
    input.atoms = (coordgen_atom_span_t){atoms, 2, 0};
    input.bonds = (coordgen_bond_span_t){&bond, 1, 0};
    if (coordgen_probe_generate(&input, &probe) != COORDGEN_OK) return 20;
    if (probe.template_mapping_count != 0 || probe.fragment_count == 0) {
        coordgen_probe_result_free(&probe);
        return 21;
    }
    for (index = 0; index < probe.fragment_count; ++index) {
        if (probe.fragments[index].template_match != COORDGEN_INVALID_INDEX ||
            probe.fragments[index].template_mapping_count != 0) {
            coordgen_probe_result_free(&probe);
            return 22;
        }
    }
    coordgen_probe_result_free(&probe);
    return probe.owner == 0 ? 0 : 23;
}

static int fragment_parent_probe_smoke(void) {
    coordgen_atom_input_t atoms[5] = {0};
    coordgen_bond_input_t bonds[4] = {0};
    coordgen_input_t input = {0};
    coordgen_probe_result_t probe = {0};
    uint32_t index;
    static const uint32_t kinds[13] = {0, 1, 2, 3, 3, 0, 1, 2, 0, 1, 2, 3, 3};
    static const uint32_t fragments[13] = {0, 0, 0, 0, 0, 1, 1, 1, 2, 2, 2, 2, 2};
    static const uint32_t counts[13] = {1, 7, 1, 2, 2, 2, 7, 5, 2, 7, 5, 2, 2};
    static const uint32_t tiers[13] = {0, 2, 3, 4, 4, 0, 2, 3, 0, 2, 3, 4, 4};
    static const uint32_t affected[13] = {0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1};
    static const uint32_t affected_starts[13] = {0, 0, 0, 0, 1, 2, 2, 2, 2, 2, 2, 2, 3};
    static const uint32_t affected_atoms[4] = {0, 1, 4, 3};
    static const uint32_t pivots[13] = {
        COORDGEN_INVALID_INDEX, COORDGEN_INVALID_INDEX, COORDGEN_INVALID_INDEX, 1, 0,
        COORDGEN_INVALID_INDEX, COORDGEN_INVALID_INDEX, COORDGEN_INVALID_INDEX,
        COORDGEN_INVALID_INDEX, COORDGEN_INVALID_INDEX, COORDGEN_INVALID_INDEX, 3, 4
    };

    for (index = 0; index < 5; ++index) {
        atoms[index].atomic_number = 6;
        atoms[index].stereo_looking_from = COORDGEN_INVALID_INDEX;
        atoms[index].stereo_atom_a = COORDGEN_INVALID_INDEX;
        atoms[index].stereo_atom_b = COORDGEN_INVALID_INDEX;
    }
    for (index = 0; index < 4; ++index) {
        bonds[index].start = index;
        bonds[index].end = index + 1;
        bonds[index].order = COORDGEN_BOND_SINGLE;
        bonds[index].stereo_atom_a = COORDGEN_INVALID_INDEX;
        bonds[index].stereo_atom_b = COORDGEN_INVALID_INDEX;
        bonds[index].crossing_penalty_multiplier = 1.0f;
    }
    input.options = coordgen_default_options();
    input.atoms = (coordgen_atom_span_t){atoms, 5, 0};
    input.bonds = (coordgen_bond_span_t){bonds, 4, 0};
    if (coordgen_probe_generate(&input, &probe) != COORDGEN_OK) return 24;
    if (probe.fragment_count != 3 ||
        probe.fragments[0].parent != COORDGEN_INVALID_INDEX ||
        probe.fragments[1].parent != 0 ||
        probe.fragments[2].parent != 1) {
        coordgen_probe_result_free(&probe);
        return 25;
    }
    if (probe.dof_count != 13 || probe.dof_affected_atoms.len != 4) {
        coordgen_probe_result_free(&probe);
        return 26;
    }
    for (index = 0; index < 4; ++index) {
        if (probe.dof_affected_atoms.ptr[index] != affected_atoms[index]) {
            coordgen_probe_result_free(&probe);
            return 27;
        }
    }
    for (index = 0; index < probe.dof_count; ++index) {
        const coordgen_probe_dof_t *d = &probe.dofs[index];
        const float expected_penalty = (index == 5 || index == 8) ? 10.0f : 0.0f;
        if (d->id != index || d->kind != kinds[index] || d->fragment != fragments[index] ||
            d->current_state != 0 || d->optimal_state != 0 || d->state_count != counts[index] ||
            d->tier != tiers[index] || d->affected_start != affected_starts[index] ||
            d->affected_count != affected[index] ||
            d->atom_a != pivots[index] || d->atom_b != COORDGEN_INVALID_INDEX ||
            d->current_penalty != expected_penalty) {
            coordgen_probe_result_free(&probe);
            return 28;
        }
    }
    coordgen_probe_result_free(&probe);
    return probe.owner == 0 ? 0 : 29;
}

/* cgz-7v2.8 regression: build_from_fragments must not be silently ignored.
 * Upstream's sketcherMinimizer::buildFromFragments(bool) is an imperative
 * pipeline step (forwards to CoordgenMinimizer::buildFromFragments(bool
 * firstTime) const), not a stored option; runGenerateCoordinates() already
 * calls it unconditionally with firstTime=true, and no upstream call site
 * ever passes false. There is no faithful way for this ABI flag to change
 * coordgen_generate()'s single-shot pipeline, so a nonzero value must be
 * rejected rather than silently accepted-and-ignored. Before the fix, the
 * adapter called buildFromFragments() before minimizer.initialize() (a
 * guaranteed no-op, since m_molecules is empty at that point) and returned
 * COORDGEN_OK regardless of the flag; this check fails against that code
 * because it observes COORDGEN_OK instead of COORDGEN_ERROR_UNSUPPORTED. */
static int build_from_fragments_rejected_smoke(void) {
    coordgen_atom_input_t atoms[3] = {0};
    coordgen_bond_input_t bonds[3] = {0};
    coordgen_input_t input = {0};
    coordgen_result_t result = {0};
    coordgen_error_t error;
    uint32_t index;

    for (index = 0; index < 3; ++index) {
        atoms[index].atomic_number = 6;
        atoms[index].stereo_looking_from = COORDGEN_INVALID_INDEX;
        atoms[index].stereo_atom_a = COORDGEN_INVALID_INDEX;
        atoms[index].stereo_atom_b = COORDGEN_INVALID_INDEX;
    }
    bonds[0] = (coordgen_bond_input_t){0, 1, COORDGEN_BOND_SINGLE, 0, 0,
                                        COORDGEN_INVALID_INDEX, COORDGEN_INVALID_INDEX, 0, 1.0f, 0};
    bonds[1] = (coordgen_bond_input_t){1, 2, COORDGEN_BOND_SINGLE, 0, 0,
                                        COORDGEN_INVALID_INDEX, COORDGEN_INVALID_INDEX, 0, 1.0f, 0};
    bonds[2] = (coordgen_bond_input_t){2, 0, COORDGEN_BOND_SINGLE, 0, 0,
                                        COORDGEN_INVALID_INDEX, COORDGEN_INVALID_INDEX, 0, 1.0f, 0};
    input.options = coordgen_default_options();
    input.options.build_from_fragments = 1;
    input.atoms = (coordgen_atom_span_t){atoms, 3, 0};
    input.bonds = (coordgen_bond_span_t){bonds, 3, 0};

    error = coordgen_generate(&input, &result);
    if (error == COORDGEN_OK) {
        coordgen_result_free(&result);
        return 30;
    }
    if (error != COORDGEN_ERROR_UNSUPPORTED) return 31;
    if (result.owner != 0) return 32;

    /* The default (0/false) must remain accepted: it matches actual pinned
     * behavior, since fragments are always built during generation. */
    input.options.build_from_fragments = 0;
    if (coordgen_generate(&input, &result) != COORDGEN_OK) return 33;
    coordgen_result_free(&result);
    return 0;
}

int main(void) {
    const int stable = stable_api_smoke();
    const int probe = stable != 0 ? stable : probe_api_smoke();
    const int no_template = probe != 0 ? probe : no_template_probe_smoke();
    const int fragment_parent = no_template != 0 ? no_template : fragment_parent_probe_smoke();
    return fragment_parent != 0 ? fragment_parent : build_from_fragments_rejected_smoke();
}
