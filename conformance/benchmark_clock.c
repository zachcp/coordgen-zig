#define _POSIX_C_SOURCE 200809L

#include <stdint.h>
#include <time.h>

uint64_t cgz_benchmark_now_ns(void) {
    struct timespec now;
    if (clock_gettime(CLOCK_MONOTONIC, &now) != 0) {
        return 0;
    }
    return (uint64_t)now.tv_sec * UINT64_C(1000000000) + (uint64_t)now.tv_nsec;
}
