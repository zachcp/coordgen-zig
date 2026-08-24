#include <cstdint>
#include <cstring>
#include <new>

#ifndef CGZ_ALLOCATOR_DESCENDING
#define CGZ_ALLOCATOR_DESCENDING 1
#endif

int main()
{
    void* first = ::operator new(0);
    void* second = ::operator new(0);
    if (first == nullptr || second == nullptr || first == second) return 1;
#if CGZ_ALLOCATOR_DESCENDING
    if (reinterpret_cast<std::uintptr_t>(first) <= reinterpret_cast<std::uintptr_t>(second)) return 3;
#else
    if (reinterpret_cast<std::uintptr_t>(first) >= reinterpret_cast<std::uintptr_t>(second)) return 3;
#endif
    std::memset(first, 0xA5, 1);
    std::memset(second, 0x5A, 1);
    ::operator delete(first);
    ::operator delete(second);

    constexpr std::size_t alignment = 256;
    void* aligned = ::operator new(1, std::align_val_t(alignment));
    if (aligned == nullptr ||
        reinterpret_cast<std::uintptr_t>(aligned) % alignment != 0) {
        return 2;
    }
    std::memset(aligned, 0x3C, 1);
    ::operator delete(aligned, std::align_val_t(alignment));
    return 0;
}
