/*
 * Out-of-tree C++ consumer for the installed C ABI. This deliberately uses
 * only the installed C header and archive; the optional C++ facade does not
 * exist yet, so C++ interoperability means the stable extern "C" surface is
 * usable with ordinary C++17 types and initialization.
 */
#include "coordgen_abi.h"

#include <type_traits>

static_assert(std::is_standard_layout<coordgen_input_t>::value,
              "input DTO must remain standard-layout in C++");
static_assert(std::is_standard_layout<coordgen_result_t>::value,
              "result DTO must remain standard-layout in C++");

int main() {
    coordgen_input_t input = {};
    coordgen_result_t result = {};
    input.options = coordgen_default_options();

    if (input.options.precision != COORDGEN_PRECISION_STANDARD) return 1;
    if (coordgen_generate(&input, &result) != COORDGEN_ERROR_EMPTY_GRAPH) return 2;
    if (result.owner != nullptr || result.coordinates.ptr != nullptr) return 3;

    coordgen_result_free(&result);
    return result.owner == nullptr && result.coordinates.ptr == nullptr ? 0 : 4;
}
