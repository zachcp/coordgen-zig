#ifndef COORDGEN_PROBE_H
#define COORDGEN_PROBE_H

#include <stdint.h>

/* Conformance-only. Never install this header or promise version stability. */
#define COORDGEN_PROBE_VERSION UINT32_C(1)

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

#endif /* COORDGEN_PROBE_H */
