#include "coordgen_probe.h"

#include <stddef.h>

_Static_assert(sizeof(coordgen_probe_dof_t) == 56, "DOF probe layout");
_Static_assert(_Alignof(coordgen_probe_dof_t) == 4, "DOF probe alignment");
_Static_assert(offsetof(coordgen_probe_dof_t, current_penalty) == 48,
               "DOF probe penalty offset");
_Static_assert(COORDGEN_PROBE_DOF_FLIP_RING == 6, "DOF probe variants");
_Static_assert(sizeof(coordgen_probe_ring_t) == 8, "ring probe layout");
_Static_assert(sizeof(coordgen_probe_template_mapping_t) == 8,
               "template mapping keeps both input and template atom IDs");
_Static_assert(sizeof(coordgen_probe_fragment_t) == 48, "fragment probe layout");
_Static_assert(sizeof(coordgen_probe_component_t) == 40, "component probe layout");
_Static_assert(offsetof(coordgen_probe_fragment_t, template_match) == 36,
               "template match remains explicit in the unstable probe");
_Static_assert(offsetof(coordgen_probe_component_t, transform) == 16,
               "component transform offset");
