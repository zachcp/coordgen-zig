#define _POSIX_C_SOURCE 200809L

#include <stdint.h>
#include <time.h>

static uint64_t clock_now_ns(clockid_t clock_id) {
    struct timespec now;
    if (clock_gettime(clock_id, &now) != 0) {
        return 0;
    }
    return (uint64_t)now.tv_sec * UINT64_C(1000000000) + (uint64_t)now.tv_nsec;
}

uint64_t cgz_benchmark_now_ns(void) {
    return clock_now_ns(CLOCK_MONOTONIC);
}

uint64_t cgz_benchmark_thread_cpu_ns(void) {
    return clock_now_ns(CLOCK_THREAD_CPUTIME_ID);
}
