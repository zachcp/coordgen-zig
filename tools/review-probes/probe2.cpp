// Same probe, but with a global operator new that hands out addresses in
// DESCENDING order, inverting every pointer comparison used by std::set /
// std::map keyed on object pointers.
#include <cstdio>
#include <cstdlib>
#include <cstddef>
#include <new>
#include <string>
#include <sys/mman.h>
#include <vector>
#include "sketcherMinimizer.h"
#include "test/coordgenBasicSMILES.h"

namespace
{
constexpr size_t POOL = 1ull << 30; // 1 GiB
char* pool_base = nullptr;
char* pool_cur = nullptr;
bool descending = false;

void init_pool()
{
    void* p = mmap(nullptr, POOL, PROT_READ | PROT_WRITE,
                   MAP_PRIVATE | MAP_ANON, -1, 0);
    if (p == MAP_FAILED) {
        abort();
    }
    pool_base = static_cast<char*>(p);
    pool_cur = pool_base + POOL;
}
} // namespace

void* operator new(size_t n)
{
    if (!descending) {
        void* p = malloc(n ? n : 1);
        if (!p) throw std::bad_alloc();
        return p;
    }
    if (!pool_base) init_pool();
    size_t sz = (n + 63) & ~size_t(63);
    pool_cur -= sz;
    if (pool_cur < pool_base) abort();
    return pool_cur;
}

void operator delete(void* p) noexcept
{
    if (pool_base && p >= pool_base && p < pool_base + POOL) return; // leak
    free(p);
}
void operator delete(void* p, size_t) noexcept { operator delete(p); }
void* operator new[](size_t n) { return operator new(n); }
void operator delete[](void* p) noexcept { operator delete(p); }
void operator delete[](void* p, size_t) noexcept { operator delete(p); }

using namespace schrodinger;

static const char* SMILES[] = {
    "C1CC2CCC3CCCC4CCC(C1)C2C34",
    "CC12CCC3C(CCC4CC(O)CCC34C)C1CCC2O",
    "C1CCCCCCCCCCCCCCCC1",
    "C1CCCCC2CCCCC(CCCCC1)C2",
    "CCC(CC)C(=O)NC(CS)C(=O)NCCCN",
    "C1CC2CCC1CC2",
    "C1CCC(CC1)C1CCC(CC1)C1CCC(CC1)C1CCCCC1",
};

int main(int argc, char** argv)
{
    descending = (argc > 1 && argv[1][0] == 'd');
    const sketcherMinimizerAtom* first = nullptr;
    const sketcherMinimizerAtom* second = nullptr;
    for (auto s : SMILES) {
        auto* mol = approxSmilesParse(s);
        sketcherMinimizer minimizer;
        std::vector<sketcherMinimizerAtom*> ats = mol->getAtoms();
        if (!first && ats.size() > 1) {
            first = ats[0];
            second = ats[1];
        }
        minimizer.initialize(mol);
        minimizer.runGenerateCoordinates();
        printf("%s\n", s);
        for (size_t i = 0; i < ats.size(); ++i) {
            union { float f; unsigned int u; } x, y;
            x.f = ats[i]->getCoordinates().x();
            y.f = ats[i]->getCoordinates().y();
            printf("  %3zu %08x %08x  %.9g %.9g\n", i, x.u, y.u, x.f, y.f);
        }
    }
    fprintf(stderr, "mode=%s atom0<atom1 = %d\n", descending ? "descending" : "normal",
            (int) (first < second));
    return 0;
}
