#ifndef COORDGEN_PROBE_H
#define COORDGEN_PROBE_H

#include "coordgen_abi.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Conformance-only. Never install this header or promise version stability. */
#define COORDGEN_PROBE_VERSION UINT32_C(2)

typedef uint32_t coordgen_probe_dof_kind_t;
enum {
    COORDGEN_PROBE_DOF_FLIP_FRAGMENT = 0,
    COORDGEN_PROBE_DOF_CHANGE_PARENT_BOND_LENGTH = 1,
    COORDGEN_PROBE_DOF_ROTATE_FRAGMENT = 2,
    COORDGEN_PROBE_DOF_SCALE_ATOMS = 3,
    COORDGEN_PROBE_DOF_SCALE_FRAGMENT = 4,
    COORDGEN_PROBE_DOF_INVERT_BOND = 5,
    COORDGEN_PROBE_DOF_FLIP_RING = 6
};

typedef struct coordgen_probe_dof {
    uint32_t id;
    coordgen_probe_dof_kind_t kind;
    uint32_t fragment;
    uint32_t current_state;
    uint32_t optimal_state;
    uint32_t state_count;
    uint32_t tier;
    uint32_t affected_start;
    uint32_t affected_count;
    uint32_t atom_a;
    uint32_t atom_b;
    uint32_t ring;
    float current_penalty;
    int32_t variant_penalty_multiplier;
} coordgen_probe_dof_t;

/* Every pointer below is owned by coordgen_probe_result_t::owner and stays
 * valid until coordgen_probe_result_free().  Indices always address the input
 * atom array, except parent/fragment/component indices which address their
 * corresponding probe arrays. Span offsets address the matching named span. */
typedef struct coordgen_probe_ring {
    uint32_t atom_start;
    uint32_t atom_count;
} coordgen_probe_ring_t;

typedef struct coordgen_probe_template_mapping {
    uint32_t input_atom;
    uint32_t template_atom;
} coordgen_probe_template_mapping_t;

enum {
    COORDGEN_PROBE_FRAGMENT_FIXED = UINT32_C(1) << 0,
    COORDGEN_PROBE_FRAGMENT_TEMPLATED = UINT32_C(1) << 1,
    COORDGEN_PROBE_FRAGMENT_CONSTRAINED = UINT32_C(1) << 2,
    COORDGEN_PROBE_FRAGMENT_CONSTRAINED_FLIP = UINT32_C(1) << 3,
    COORDGEN_PROBE_FRAGMENT_CHAIN = UINT32_C(1) << 4
};

typedef struct coordgen_probe_fragment {
    uint32_t parent;
    uint32_t component;
    uint32_t atom_start;
    uint32_t atom_count;
    uint32_t ring_start;
    uint32_t ring_count;
    uint32_t dof_start;
    uint32_t dof_count;
    uint32_t flags;
    /* The exact embedded-template index and entries in result.template_mapping.
     * Use COORDGEN_INVALID_INDEX and an empty range for no template match. */
    uint32_t template_match;
    uint32_t template_mapping_start;
    uint32_t template_mapping_count;
} coordgen_probe_fragment_t;

typedef uint32_t coordgen_probe_transform_status_t;
enum {
    /* The oracle hook did not observe a transform for this component. */
    COORDGEN_PROBE_TRANSFORM_UNOBSERVED = 0,
    /* Captured around upstream's component placement phase. */
    COORDGEN_PROBE_TRANSFORM_OBSERVED = 1
};

typedef struct coordgen_probe_component {
    uint32_t atom_start;
    uint32_t atom_count;
    coordgen_probe_transform_status_t transform_status;
    uint32_t reserved;
    /* [m00, m01, m10, m11, tx, ty]; zero when status is UNOBSERVED. */
    float transform[6];
} coordgen_probe_component_t;

typedef struct coordgen_probe_result {
    coordgen_u32_span_t input_to_internal;
    coordgen_u32_span_t internal_to_input;
    coordgen_u32_span_t morgan_ranks;
    coordgen_u32_span_t ring_atoms;
    coordgen_u32_span_t fragment_atoms;
    coordgen_u32_span_t fragment_rings;
    coordgen_u32_span_t component_atoms;
    coordgen_u32_span_t dof_affected_atoms;
    coordgen_probe_template_mapping_t *template_mapping;
    uint32_t template_mapping_count;
    uint32_t template_mapping_reserved;
    coordgen_probe_ring_t *rings;
    uint32_t ring_count;
    uint32_t ring_reserved;
    coordgen_probe_fragment_t *fragments;
    uint32_t fragment_count;
    uint32_t fragment_reserved;
    coordgen_probe_dof_t *dofs;
    uint32_t dof_count;
    uint32_t dof_reserved;
    coordgen_probe_component_t *components;
    uint32_t component_count;
    uint32_t clean_pose;
    uint32_t reserved;
    void *owner;
} coordgen_probe_result_t;

/* Uses the same borrowed input contract as coordgen_generate. The probe is
 * conformance-only: it is built only with -Denable-oracle=true and is never
 * installed or included by the production package. */
coordgen_error_t coordgen_probe_generate(const coordgen_input_t *input,
                                         coordgen_probe_result_t *result);
void coordgen_probe_result_free(coordgen_probe_result_t *result);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* COORDGEN_PROBE_H */
