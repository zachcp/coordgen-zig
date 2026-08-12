/*
 * Compiled and run by tools/check-install-isolation against a *scratch*
 * install prefix (`-I prefix/include -L prefix/lib -lcoordgen`), never
 * against this repository's build tree or source paths. Its only job is to
 * prove the installed artifacts are self-sufficient and behave as
 * documented for a consumer who has nothing but `$prefix`.
 */
#include "coordgen_abi.h"

#include <stddef.h>

int main(void) {
    coordgen_input_t input;
    coordgen_result_t result;

    input.options = coordgen_default_options();
    input.atoms.ptr = NULL;
    input.atoms.len = 0;
    input.atoms.reserved = 0;
    input.bonds.ptr = NULL;
    input.bonds.len = 0;
    input.bonds.reserved = 0;
    input.residues.ptr = NULL;
    input.residues.len = 0;
    input.residues.reserved = 0;
    input.residue_interactions.ptr = NULL;
    input.residue_interactions.len = 0;
    input.residue_interactions.reserved = 0;
    input.extra_bonds.ptr = NULL;
    input.extra_bonds.len = 0;
    input.extra_bonds.reserved = 0;

    if (coordgen_generate(&input, &result) != COORDGEN_ERROR_EMPTY_GRAPH) {
        return 1;
    }
    if (result.owner != NULL) {
        return 1;
    }
    /* Documented safe on a zeroed failure result; also proves
     * coordgen_result_free itself resolved from the installed archive. */
    coordgen_result_free(&result);

    return 0;
}
