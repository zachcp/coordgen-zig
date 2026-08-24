/* Deterministic heap-address-order control for the oracle stability
 * classification.
 *
 * Upstream keys 94 std::set/std::map declarations on pointers to
 * sketcherMinimizerAtom, Ring, Molecule, and Fragment, so their iteration
 * order is heap address order. The baseline and architecture runners define
 * CGZ_ALLOCATOR_DESCENDING=0 and receive ascending addresses. The adversarial
 * runner defines it to 1 and receives descending addresses, inverting every
 * one of those comparisons without changing upstream source: `atom0 < atom1`
 * becomes false where it was true.
 *
 * Both sides must use this translation unit. A platform allocator does not
 * promise ascending addresses, so comparing it with the descending allocator
 * made the measured member set depend on the host and residual heap layout.
 *
 * Nothing is ever returned to the pool, because reuse would break monotonic
 * address order. The reservation is therefore sized for a whole corpus run:
 * address space, not memory. A
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

#ifndef CGZ_ALLOCATOR_DESCENDING
#define CGZ_ALLOCATOR_DESCENDING 1
#endif

namespace
{

/* One large reservation, consumed in the selected direction. The bump pointer
 * never wraps or reuses, so the reservation has to cover every allocation a
 * corpus run makes. */
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
#if CGZ_ALLOCATOR_DESCENDING
    pool_cursor = pool_base + kPoolBytes;
#else
    pool_cursor = pool_base;
#endif
}

void* allocateMonotonic(std::size_t size, std::size_t alignment)
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
#if CGZ_ALLOCATOR_DESCENDING
    if (cursor - base < requested) allocationFailure("descending pool exhausted");
    const std::uintptr_t next = (cursor - requested) & ~(alignment - 1);
    if (next < base) allocationFailure("descending pool exhausted after alignment");
    pool_cursor = reinterpret_cast<char*>(next);
    return pool_cursor;
#else
    const std::uintptr_t aligned = (cursor + alignment - 1) & ~(alignment - 1);
    const std::uintptr_t end = base + kPoolBytes;
    if (aligned < cursor || aligned > end || requested > end - aligned) {
        allocationFailure("ascending pool exhausted");
    }
    pool_cursor = reinterpret_cast<char*>(aligned + requested);
    return reinterpret_cast<void*>(aligned);
#endif
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
    return allocateMonotonic(size, alignof(std::max_align_t));
}

void* operator new(std::size_t size, std::align_val_t alignment)
{
    return allocateMonotonic(size, static_cast<std::size_t>(alignment));
}

void operator delete(void* pointer) noexcept
{
    /* Pool memory is never reused: reuse would hand out an address that is no
     * longer monotonic. The pool is released with the process. */
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
