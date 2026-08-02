/* Dumps the drug-like corpus partition as upstream's own SMILES parser builds
 * it.
 *
 * The drug-like structures are committed as tables in
 * src/conformance/corpus.zig so that the corpus is one self-contained Zig
 * artifact and does not depend on a C++ parser at generation time. This
 * program exists so that the committed tables are checked rather than
 * trusted: the build runs it against the pinned oracle and diffs its output
 * against the Zig generator's dump of the same partition.
 *
 * Conformance-only. Never linked into an installed artifact.
 */

#include <cstdio>
#include <string>
#include <vector>

#include "sketcherMinimizerAtom.h"
#include "sketcherMinimizerBond.h"
#include "sketcherMinimizerMolecule.h"
#include "test/coordgenBasicSMILES.h"

namespace
{

/* Pinned drug-like corpus. Taken from the review probe behind cgz-r05, whose
 * measurements established that realistic structures are bit-identical across
 * architectures and optimization levels: two fused polycyclics, a steroid
 * skeleton, a seventeen-membered macrocycle, a bridged bicyclic pair, a
 * peptide-like chain, and a four-ring biaryl chain. */
const char* const kSmiles[] = {
    "C1CC2CCC3CCCC4CCC(C1)C2C34",
    "CC12CCC3C(CCC4CC(O)CCC34C)C1CCC2O",
    "C1CCCCCCCCCCCCCCCC1",
    "C1CCCCC2CCCCC(CCCCC1)C2",
    "CCC(CC)C(=O)NC(CS)C(=O)NCCCN",
    "C1CC2CCC1CC2",
    "C1CCC(CC1)C1CCC(CC1)C1CCC(CC1)C1CCCCC1",
};

} // namespace

int main()
{
    unsigned index = 0;
    for (const char* smiles : kSmiles) {
        auto* molecule = schrodinger::approxSmilesParse(std::string(smiles));
        const std::vector<sketcherMinimizerAtom*>& atoms = molecule->getAtoms();
        const std::vector<sketcherMinimizerBond*>& bonds = molecule->getBonds();

        std::printf("molecule drug_like %u atoms=%zu bonds=%zu\n", index,
                    atoms.size(), bonds.size());
        for (std::size_t atom = 0; atom < atoms.size(); ++atom) {
            std::printf("atom %zu z=%d q=%d\n", atom,
                        atoms[atom]->getAtomicNumber(), atoms[atom]->charge);
        }
        for (std::size_t bond = 0; bond < bonds.size(); ++bond) {
            std::size_t start = 0;
            std::size_t end = 0;
            for (std::size_t atom = 0; atom < atoms.size(); ++atom) {
                if (atoms[atom] == bonds[bond]->getStartAtom()) start = atom;
                if (atoms[atom] == bonds[bond]->getEndAtom()) end = atom;
            }
            std::printf("bond %zu %zu %zu order=%d\n", bond, start, end,
                        bonds[bond]->getBondOrder());
        }
        ++index;
    }
    return 0;
}
