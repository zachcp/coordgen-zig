#ifndef COORDGEN_ORACLE_HOOK_HPP
#define COORDGEN_ORACLE_HOOK_HPP

#include <cstddef>
#include <vector>

class sketcherMinimizerAtom;
class sketcherMinimizerFragment;
class sketcherMinimizerMolecule;

/* This header is included only by build-generated, patched copies of two
 * pinned oracle translation units.  The installed native library never sees
 * it, and the hook functions are implemented by oracle_adapter.cpp. */
namespace coordgen_oracle_hook {

class CaptureSink {
  public:
    virtual ~CaptureSink() = default;

    virtual void recordTemplateMatch(
        sketcherMinimizerFragment* fragment,
        const std::vector<sketcherMinimizerAtom*>& atoms,
        std::size_t template_index,
        const std::vector<unsigned int>& mapping) = 0;

    virtual void captureComponentTransforms(
        const std::vector<sketcherMinimizerMolecule*>& components,
        bool after_placement) = 0;
};

void setCaptureSink(CaptureSink* sink);
void recordTemplateMatch(sketcherMinimizerFragment* fragment,
                         const std::vector<sketcherMinimizerAtom*>& atoms,
                         std::size_t template_index,
                         const std::vector<unsigned int>& mapping);

void captureComponentTransforms(
    const std::vector<sketcherMinimizerMolecule*>& components,
    bool after_placement);

} // namespace coordgen_oracle_hook

#endif /* COORDGEN_ORACLE_HOOK_HPP */
