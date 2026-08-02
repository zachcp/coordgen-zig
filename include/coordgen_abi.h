#ifndef COORDGEN_ABI_H
#define COORDGEN_ABI_H

#include <stdint.h>
#include <string.h>

#ifdef __cplusplus
extern "C" {
#endif

#define COORDGEN_ABI_VERSION UINT32_C(1)
#define COORDGEN_INVALID_INDEX UINT32_MAX
#define COORDGEN_BOND_LENGTH 50.0f
#define COORDGEN_PRECISION_QUICK 0.2f
#define COORDGEN_PRECISION_STANDARD 1.0f
#define COORDGEN_PRECISION_BEST 3.0f

typedef uint32_t coordgen_error_t;
enum {
    COORDGEN_OK = 0,
    COORDGEN_ERROR_EMPTY_GRAPH = 1,
    COORDGEN_ERROR_TOO_MANY_ITEMS = 2,
    COORDGEN_ERROR_INVALID_ATOMIC_NUMBER = 3,
    COORDGEN_ERROR_INVALID_BOND_ORDER = 4,
    COORDGEN_ERROR_INVALID_ATOM_INDEX = 5,
    COORDGEN_ERROR_INVALID_STEREO = 6,
    COORDGEN_ERROR_INVALID_COORDINATE = 7,
    COORDGEN_ERROR_INVALID_OPTION = 8,
    COORDGEN_ERROR_INVALID_MAPPING = 9,
    COORDGEN_ERROR_OUT_OF_MEMORY = 10,
    COORDGEN_ERROR_UNSUPPORTED = 11,
    COORDGEN_ERROR_INTERNAL = 12
};

typedef uint32_t coordgen_bond_order_t;
enum {
    COORDGEN_BOND_ZERO = 0,
    COORDGEN_BOND_SINGLE = 1,
    COORDGEN_BOND_DOUBLE = 2,
    COORDGEN_BOND_TRIPLE = 3
};

typedef uint32_t coordgen_atom_stereo_t;
enum {
    COORDGEN_ATOM_STEREO_UNSPECIFIED = 0,
    COORDGEN_ATOM_STEREO_CLOCKWISE = 1,
    COORDGEN_ATOM_STEREO_COUNTER_CLOCKWISE = 2,
    COORDGEN_ATOM_STEREO_R = 3,
    COORDGEN_ATOM_STEREO_S = 4
};

typedef uint32_t coordgen_bond_stereo_t;
enum {
    COORDGEN_BOND_STEREO_UNSPECIFIED = 0,
    COORDGEN_BOND_STEREO_CIS = 1,
    COORDGEN_BOND_STEREO_TRANS = 2,
    COORDGEN_BOND_STEREO_Z = 3,
    COORDGEN_BOND_STEREO_E = 4
};

typedef uint32_t coordgen_bond_display_t;
enum {
    COORDGEN_BOND_DISPLAY_NONE = 0,
    COORDGEN_BOND_DISPLAY_SOLID_FORWARD = 1,
    COORDGEN_BOND_DISPLAY_SOLID_REVERSE = 2,
    COORDGEN_BOND_DISPLAY_HASHED_FORWARD = 3,
    COORDGEN_BOND_DISPLAY_HASHED_REVERSE = 4
};

enum {
    COORDGEN_ATOM_HAS_TEMPLATE_COORDINATES = UINT32_C(1) << 0,
    COORDGEN_ATOM_HAS_COORDINATES_3D = UINT32_C(1) << 1,
    COORDGEN_ATOM_FIXED = UINT32_C(1) << 2,
    COORDGEN_ATOM_CONSTRAINED = UINT32_C(1) << 3,
    COORDGEN_ATOM_HIDDEN = UINT32_C(1) << 4,
    COORDGEN_BOND_SKIP = UINT32_C(1) << 0
};

typedef struct coordgen_vec2 {
    float x;
    float y;
} coordgen_vec2_t;

typedef struct coordgen_vec3 {
    float x;
    float y;
    float z;
} coordgen_vec3_t;

typedef struct coordgen_string_view {
    const uint8_t *ptr;
    uint32_t len;
    uint32_t reserved;
} coordgen_string_view_t;

typedef struct coordgen_index_span {
    const uint32_t *ptr;
    uint32_t len;
    uint32_t reserved;
} coordgen_index_span_t;

typedef struct coordgen_atom_input {
    uint32_t atomic_number;
    int32_t formal_charge;
    uint32_t flags;
    coordgen_atom_stereo_t stereo;
    uint32_t stereo_looking_from;
    uint32_t stereo_atom_a;
    uint32_t stereo_atom_b;
    coordgen_vec2_t template_coordinates;
    coordgen_vec3_t coordinates_3d;
    uint32_t reserved;
} coordgen_atom_input_t;

typedef struct coordgen_bond_input {
    uint32_t start;
    uint32_t end;
    coordgen_bond_order_t order;
    uint32_t flags;
    coordgen_bond_stereo_t stereo;
    uint32_t stereo_atom_a;
    uint32_t stereo_atom_b;
    coordgen_bond_display_t display;
    float crossing_penalty_multiplier;
    uint32_t reserved;
} coordgen_bond_input_t;

typedef struct coordgen_residue_input {
    uint32_t atom;
    int32_t residue_number;
    uint32_t closest_ligand_atom;
    uint32_t reserved;
    coordgen_string_view_t chain;
} coordgen_residue_input_t;

typedef struct coordgen_residue_interaction_input {
    uint32_t start;
    uint32_t end;
    coordgen_index_span_t other_start_atoms;
    coordgen_index_span_t other_end_atoms;
    float crossing_penalty_multiplier;
    uint32_t reserved[3];
} coordgen_residue_interaction_input_t;

#define COORDGEN_DECLARE_CONST_SPAN(name, element_type) \
    typedef struct name {                              \
        const element_type *ptr;                       \
        uint32_t len;                                  \
        uint32_t reserved;                             \
    } name##_t

#define COORDGEN_DECLARE_MUT_SPAN(name, element_type) \
    typedef struct name {                            \
        element_type *ptr;                           \
        uint32_t len;                                \
        uint32_t reserved;                           \
    } name##_t

COORDGEN_DECLARE_CONST_SPAN(coordgen_atom_span, coordgen_atom_input_t);
COORDGEN_DECLARE_CONST_SPAN(coordgen_bond_span, coordgen_bond_input_t);
COORDGEN_DECLARE_CONST_SPAN(coordgen_residue_span, coordgen_residue_input_t);
COORDGEN_DECLARE_CONST_SPAN(coordgen_residue_interaction_span,
                            coordgen_residue_interaction_input_t);
COORDGEN_DECLARE_MUT_SPAN(coordgen_vec2_span, coordgen_vec2_t);
COORDGEN_DECLARE_MUT_SPAN(coordgen_u32_span, uint32_t);

#undef COORDGEN_DECLARE_CONST_SPAN
#undef COORDGEN_DECLARE_MUT_SPAN

typedef struct coordgen_options {
    float precision;
    uint32_t score_residue_interactions;
    uint32_t treat_nonterminal_bonds_to_metal_as_zero_order;
    uint32_t even_angles;
    uint32_t skip_minimization;
    uint32_t force_open_macrocycles;
    uint32_t constrain_all_atoms;
    uint32_t build_from_fragments;
    uint32_t debug_coordinates;
    uint32_t load_templates;
    coordgen_string_view_t template_directory;
} coordgen_options_t;

typedef struct coordgen_input {
    coordgen_options_t options;
    coordgen_atom_span_t atoms;
    coordgen_bond_span_t bonds;
    coordgen_residue_span_t residues;
    coordgen_residue_interaction_span_t residue_interactions;
    coordgen_bond_span_t extra_bonds;
} coordgen_input_t;

typedef struct coordgen_result {
    coordgen_vec2_span_t coordinates;
    coordgen_u32_span_t input_to_internal;
    coordgen_u32_span_t internal_to_input;
    coordgen_u32_span_t effective_bond_orders;
    coordgen_u32_span_t bond_displays;
    coordgen_u32_span_t atom_stereo;
    uint32_t clean_pose;
    uint32_t reserved;
    void *owner;
} coordgen_result_t;

static inline coordgen_options_t coordgen_default_options(void) {
    coordgen_options_t value;
    memset(&value, 0, sizeof value);
    value.precision = COORDGEN_PRECISION_STANDARD;
    value.score_residue_interactions = 1;
    value.treat_nonterminal_bonds_to_metal_as_zero_order = 1;
    value.load_templates = 1;
    return value;
}

/*
 * Input pointers are borrowed only for this call. On success, result is owned
 * by the caller and must be passed exactly once to coordgen_result_free(). On
 * failure, result is zeroed and requires no cleanup. Arrays in result use
 * caller input order except the explicitly named permutation maps.
 */
coordgen_error_t coordgen_generate(const coordgen_input_t *input,
                                   coordgen_result_t *result);
void coordgen_result_free(coordgen_result_t *result);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* COORDGEN_ABI_H */
