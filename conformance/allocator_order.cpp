/* Heap-address-order adversary for the oracle stability classification.
 *
 * Upstream keys 94 std::set/std::map declarations on pointers to
 * sketcherMinimizerAtom, Ring, Molecule, and Fragment, so their iteration
 * order is heap address order. An allocator that hands out *descending* addresses
 * inverts every one of those comparisons without changing a line of upstream
 * source: `atom0 < atom1` becomes false where it was true.
 *
 * Linking this translation unit into an oracle build is therefore the whole
 * pointer-order axis of the 2x2 classification. It is linked only into the
 * descending variant; the ascending variant keeps the platform allocator.
 *
 * Nothing is ever returned to the pool, because reusing an address would hand
 * out one that is no longer monotonically descending. The reservation is
 * therefore sized for a whole corpus run: address space, not memory. A
 * 2000-member adversarial run walks about 8 GiB of addresses while touching
 * roughly 330 MiB of pages, so the reservation is large and the resident set
 * is not.
 *
 * Conformance-only. Never linked into an installed artifact.
 */

#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <new>
#include <sys/mman.h>

namespace
{

/* One large reservation, consumed downward. The bump pointer never wraps or
 * reuses - either behaviour would reintroduce ascending addresses - so the
 * reservation has to cover every allocation a corpus run makes. */
constexpr std::size_t kPoolBytes = 1ull << 35;
constexpr std::size_t kAlignment = 64;

char* pool_base = nullptr;
char* pool_cursor = nullptr;

[[noreturn]] void allocationFailure(const char* reason)
{
    std::fprintf(stderr, "allocator_order: %s\n", reason);
    std::abort();
}

void initPool()
{
    int flags = MAP_PRIVATE | MAP_ANON;
#ifdef MAP_NORESERVE
    /* The reservation is address space; only the pages actually written are
     * ever backed. */
    flags |= MAP_NORESERVE;
#endif
    void* reservation = mmap(nullptr, kPoolBytes, PROT_READ | PROT_WRITE, flags, -1, 0);
    if (reservation == MAP_FAILED) std::abort();
    pool_base = static_cast<char*>(reservation);
    pool_cursor = pool_base + kPoolBytes;
}

void* allocateDescending(std::size_t size, std::size_t alignment)
{
    if (pool_base == nullptr) initPool();
    if (alignment < kAlignment) alignment = kAlignment;
    if ((alignment & (alignment - 1)) != 0) allocationFailure("non-power-of-two alignment");
    const std::size_t requested = size == 0 ? 1 : size;
    if (requested > std::numeric_limits<std::size_t>::max() - (alignment - 1)) {
        allocationFailure("allocation size overflow");
    }
    const std::uintptr_t base = reinterpret_cast<std::uintptr_t>(pool_base);
    const std::uintptr_t cursor = reinterpret_cast<std::uintptr_t>(pool_cursor);
    if (cursor - base < requested) allocationFailure("descending pool exhausted");
    const std::uintptr_t next = (cursor - requested) & ~(alignment - 1);
    if (next < base) allocationFailure("descending pool exhausted after alignment");
    pool_cursor = reinterpret_cast<char*>(next);
    return pool_cursor;
}

bool isPoolPointer(void* pointer)
{
    if (pool_base == nullptr || pointer == nullptr) return false;
    const std::uintptr_t base = reinterpret_cast<std::uintptr_t>(pool_base);
    const std::uintptr_t address = reinterpret_cast<std::uintptr_t>(pointer);
    return address >= base && address < base + kPoolBytes;
}

} // namespace

void* operator new(std::size_t size)
{
    return allocateDescending(size, alignof(std::max_align_t));
}

void* operator new(std::size_t size, std::align_val_t alignment)
{
    return allocateDescending(size, static_cast<std::size_t>(alignment));
}

void operator delete(void* pointer) noexcept
{
    /* Pool memory is never reused: reuse would hand out an address that is no
     * longer monotonically descending. The pool is released with the process. */
    if (isPoolPointer(pointer)) return;
    std::free(pointer);
}

void operator delete(void* pointer, std::size_t) noexcept { operator delete(pointer); }
void operator delete(void* pointer, std::align_val_t) noexcept { operator delete(pointer); }
void operator delete(void* pointer, std::size_t, std::align_val_t) noexcept { operator delete(pointer); }
void* operator new[](std::size_t size) { return operator new(size); }
void* operator new[](std::size_t size, std::align_val_t alignment)
{
    return operator new(size, alignment);
}
void operator delete[](void* pointer) noexcept { operator delete(pointer); }
void operator delete[](void* pointer, std::size_t) noexcept { operator delete(pointer); }
void operator delete[](void* pointer, std::align_val_t alignment) noexcept { operator delete(pointer, alignment); }
void operator delete[](void* pointer, std::size_t size, std::align_val_t alignment) noexcept
{
    operator delete(pointer, size, alignment);
}
