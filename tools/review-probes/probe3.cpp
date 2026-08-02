// Randomized corpus + allocator-order adversary.
// Builds N random connected molecular graphs (chains, rings, fused rings,
// macrocycles, branches, heteroatoms, charges, multi-component), runs coordgen,
// and prints bit-exact coordinates. Run with "d" to invert heap address order.
#include <cstdio>
#include <cstdlib>
#include <cstddef>
#include <cstdint>
#include <new>
#include <string>
#include <sys/mman.h>
#include <vector>
#include "sketcherMinimizer.h"

namespace
{
constexpr size_t POOL = 1ull << 32;
char* pool_base = nullptr;
char* pool_cur = nullptr;
bool descending = false;

void init_pool()
{
    void* p = mmap(nullptr, POOL, PROT_READ | PROT_WRITE,
                   MAP_PRIVATE | MAP_ANON, -1, 0);
    if (p == MAP_FAILED) abort();
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
    if (pool_base && p >= pool_base && p < pool_base + POOL) return;
    free(p);
}
void operator delete(void* p, size_t) noexcept { operator delete(p); }
void* operator new[](size_t n) { return operator new(n); }
void operator delete[](void* p) noexcept { operator delete(p); }
void operator delete[](void* p, size_t) noexcept { operator delete(p); }

// xorshift64* — reproducible across platforms
static uint64_t rng_state = 0x243F6A8885A308D3ull;
static uint64_t rnd()
{
    uint64_t x = rng_state;
    x ^= x >> 12;
    x ^= x << 25;
    x ^= x >> 27;
    rng_state = x;
    return x * 0x2545F4914F6CDD1Dull;
}
static unsigned rnd_n(unsigned n) { return (unsigned) (rnd() % n); }

static const int ELEMENTS[] = {6, 6, 6, 6, 7, 8, 16, 9, 17, 15, 26, 30};

static sketcherMinimizerMolecule* random_mol(unsigned seed)
{
    rng_state = 0x9E3779B97F4A7C15ull ^ (seed * 0xBF58476D1CE4E5B9ull);
    for (int i = 0; i < 8; ++i) rnd();

    auto* mol = new sketcherMinimizerMolecule();
    std::vector<sketcherMinimizerAtom*> atoms;
    std::vector<sketcherMinimizerBond*> bonds;

    unsigned n_frag = 1 + rnd_n(3); // disconnected components
    for (unsigned f = 0; f < n_frag; ++f) {
        unsigned base = (unsigned) atoms.size();
        unsigned n = 3 + rnd_n(28);
        for (unsigned i = 0; i < n; ++i) {
            auto* a = mol->addNewAtom();
            a->molecule = mol;
            a->setAtomicNumber(ELEMENTS[rnd_n(12)]);
            if (rnd_n(12) == 0) a->charge = (int) rnd_n(3) - 1;
            atoms.push_back(a);
        }
        // spanning tree
        for (unsigned i = 1; i < n; ++i) {
            unsigned parent = base + rnd_n(i);
            auto* b = mol->addNewBond(atoms[parent], atoms[base + i]);
            unsigned r = rnd_n(10);
            b->setBondOrder(r < 6 ? 1 : (r < 9 ? 2 : (r == 9 ? 3 : 0)));
            bonds.push_back(b);
        }
        // extra ring-closing bonds (some short, some long -> macrocycles)
        unsigned extra = rnd_n(4);
        for (unsigned k = 0; k < extra && n > 3; ++k) {
            unsigned i = base + rnd_n(n);
            unsigned j = base + rnd_n(n);
            if (i == j) continue;
            bool dup = false;
            for (auto* b : bonds)
                if ((b->getStartAtom() == atoms[i] && b->getEndAtom() == atoms[j]) ||
                    (b->getStartAtom() == atoms[j] && b->getEndAtom() == atoms[i]))
                    dup = true;
            if (dup) continue;
            auto* b = mol->addNewBond(atoms[i], atoms[j]);
            b->setBondOrder(1);
            bonds.push_back(b);
        }
    }
    sketcherMinimizerMolecule::assignBondsAndNeighbors(mol->getAtoms(),
                                                       mol->getBonds());
    return mol;
}

int main(int argc, char** argv)
{
    descending = (argc > 1 && argv[1][0] == 'd');
    unsigned count = (argc > 2) ? (unsigned) atoi(argv[2]) : 200;
    for (unsigned s = 0; s < count; ++s) {
        auto* mol = random_mol(s);
        std::vector<sketcherMinimizerAtom*> ats = mol->getAtoms();
        sketcherMinimizer minimizer;
        minimizer.initialize(mol);
        bool ok = minimizer.runGenerateCoordinates();
        printf("mol %u n=%zu ok=%d\n", s, ats.size(), (int) ok);
        for (size_t i = 0; i < ats.size(); ++i) {
            union { float f; unsigned int u; } x, y;
            x.f = ats[i]->getCoordinates().x();
            y.f = ats[i]->getCoordinates().y();
            printf("  %3zu %08x %08x\n", i, x.u, y.u);
        }
    }
    return 0;
}
