#include "coordgen_probe.h"

#include <stddef.h>

_Static_assert(sizeof(coordgen_probe_dof_t) == 56, "DOF probe layout");
_Static_assert(_Alignof(coordgen_probe_dof_t) == 4, "DOF probe alignment");
_Static_assert(offsetof(coordgen_probe_dof_t, current_penalty) == 48,
               "DOF probe penalty offset");
_Static_assert(COORDGEN_PROBE_DOF_FLIP_RING == 6, "DOF probe variants");
