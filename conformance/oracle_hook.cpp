#include "coordgen_oracle_hook.hpp"

namespace coordgen_oracle_hook {

thread_local CaptureSink* active_sink = nullptr;

void setCaptureSink(CaptureSink* sink) { active_sink = sink; }

void recordTemplateMatch(sketcherMinimizerFragment* fragment,
                         const std::vector<sketcherMinimizerAtom*>& atoms,
                         std::size_t template_index,
                         const std::vector<unsigned int>& mapping)
{
    if (active_sink != nullptr) {
        active_sink->recordTemplateMatch(fragment, atoms, template_index, mapping);
    }
}

void captureComponentTransforms(
    const std::vector<sketcherMinimizerMolecule*>& components,
    bool after_placement)
{
    if (active_sink != nullptr) {
        active_sink->captureComponentTransforms(components, after_placement);
    }
}

} // namespace coordgen_oracle_hook
