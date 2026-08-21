//! Out-of-memory coverage for the port's fallible paths (cgz-5q6).
//!
//! `checkAllocationFailures` replaces `std.testing.checkAllAllocationFailures`
//! at every call site in this tree. It asserts that every allocation failure
//! propagates out, that no path leaks on the way, and that every measured run
//! reaches the allocation index it is meant to exercise.
//!
//! ## What std also asserts, and why it does not hold here
//!
//! std counts the allocations of one unconstrained run, then runs once per
//! index injecting a failure there. If a run completes without its injected
//! failure ever being requested, std reports
//! `error.NondeterministicMemoryUsage`: that run made fewer allocations than
//! the first one, so the count is not stable.
//!
//! The count is not stable here, and cannot be made so cheaply.
//! `std.ArrayList.append` grows through `Allocator.remap` and only allocates
//! when the backing allocator declines to resize in place, which is a property
//! of the heap at that moment rather than of the input. Measured on
//! `layout.basic`'s macrocycle dispatch test: the unconstrained run allocates
//! 69 and the runs at fail_index 177 and 178 of the sibling native-minimal
//! test allocate 177. Ten modules in this tree use growable containers on
//! covered paths, several with combinatorial counts where a counting pass
//! means running the enumeration twice.
//!
//! ## Making the count stable
//!
//! The allocator presented to the tested function deliberately declines every
//! in-place `resize` and `remap`. Growable containers therefore take their
//! specified allocate-copy-free fallback instead of depending on whether the
//! backing heap happens to extend a block in place. A count that still varies
//! is a property of the tested function, not of heap placement, and remains a
//! hard failure: silently skipping that index would leave later allocation
//! sites untested.
//!
//! The tests at the bottom of this file plant each violation the helper must
//! still catch: a leak on an OOM path, a swallowed OOM, and a propagated
//! non-OOM error. They are the evidence that the two kept assertions are real,
//! in the shape `tools/check-gate-strength` and its siblings use for the
//! build-step gates.

const std = @import("std");

/// Preserve allocation and free while declining in-place growth, whose success
/// depends on the backing heap's current shape. Allocator clients must already
/// support this; `realloc` falls back to allocating and copying.
const NoInPlaceAllocator = struct {
    backing: std.mem.Allocator,

    fn allocator(self: *NoInPlaceAllocator) std.mem.Allocator {
        return .{
            .ptr = self,
            .vtable = &.{
                .alloc = alloc,
                .resize = resize,
                .remap = remap,
                .free = free,
            },
        };
    }

    fn alloc(context: *anyopaque, len: usize, alignment: std.mem.Alignment, return_address: usize) ?[*]u8 {
        const self: *NoInPlaceAllocator = @ptrCast(@alignCast(context));
        return self.backing.rawAlloc(len, alignment, return_address);
    }

    fn resize(_: *anyopaque, _: []u8, _: std.mem.Alignment, _: usize, _: usize) bool {
        return false;
    }

    fn remap(_: *anyopaque, _: []u8, _: std.mem.Alignment, _: usize, _: usize) ?[*]u8 {
        return null;
    }

    fn free(context: *anyopaque, memory: []u8, alignment: std.mem.Alignment, return_address: usize) void {
        const self: *NoInPlaceAllocator = @ptrCast(@alignCast(context));
        self.backing.rawFree(memory, alignment, return_address);
    }
};

/// Run `test_fn` once per allocation it makes, failing that allocation, and
/// require that the failure propagates out as `error.OutOfMemory` with every
/// byte returned.
///
/// `test_fn`'s first parameter is the allocator; `extra_args` is a tuple
/// supplying the rest, exactly as `std.testing.checkAllAllocationFailures`
/// takes them.
pub fn checkAllocationFailures(
    backing_allocator: std.mem.Allocator,
    comptime test_fn: anytype,
    extra_args: anytype,
) !void {
    return checkAllocationFailuresInternal(backing_allocator, test_fn, extra_args, true);
}

fn checkAllocationFailuresInternal(
    backing_allocator: std.mem.Allocator,
    comptime test_fn: anytype,
    extra_args: anytype,
    report_leak: bool,
) !void {
    var no_in_place = NoInPlaceAllocator{ .backing = backing_allocator };
    const stable_allocator = no_in_place.allocator();
    // The unconstrained run has to succeed before failing anything is
    // meaningful, and it establishes how many indices there are to try.
    const total: usize = count: {
        var counting = std.testing.FailingAllocator.init(stable_allocator, .{});
        try @call(.auto, test_fn, .{counting.allocator()} ++ extra_args);
        break :count counting.alloc_index;
    };

    var fail_index: usize = 0;
    while (fail_index < total) : (fail_index += 1) {
        var failing = std.testing.FailingAllocator.init(stable_allocator, .{
            .fail_index = fail_index,
        });
        if (@call(.auto, test_fn, .{failing.allocator()} ++ extra_args)) |_| {
            // Succeeding after being denied an allocation means the failure was
            // caught and hidden. That is a real defect and still fails here.
            if (failing.has_induced_failure) return error.SwallowedOutOfMemoryError;
            // The no-in-place wrapper removes heap-layout-dependent variation.
            // A remaining count change belongs to the tested function and
            // means this allocation index was not exercised.
            return error.NondeterministicMemoryUsage;
        } else |err| switch (err) {
            error.OutOfMemory => {
                if (failing.allocated_bytes != failing.freed_bytes) {
                    if (report_leak) {
                        std.debug.print(
                            "\nleak on the out-of-memory path" ++
                                "\n  fail_index: {d} of {d}" ++
                                "\n  allocated bytes: {d}" ++
                                "\n  freed bytes: {d}" ++
                                "\n  allocations: {d}" ++
                                "\n  deallocations: {d}" ++
                                "\n  allocation that was made to fail: {f}\n",
                            .{
                                fail_index,
                                total,
                                failing.allocated_bytes,
                                failing.freed_bytes,
                                failing.allocations,
                                failing.deallocations,
                                std.debug.FormatStackTrace{
                                    .stack_trace = failing.getStackTrace(),
                                },
                            },
                        );
                    }
                    return error.MemoryLeakDetected;
                }
            },
            else => |other| return other,
        }
    }
}

fn allocatesTwiceAndCleansUp(allocator: std.mem.Allocator) !void {
    const first = try allocator.alloc(u8, 16);
    defer allocator.free(first);
    const second = try allocator.alloc(u8, 32);
    defer allocator.free(second);
}

fn leaksTheFirstAllocation(allocator: std.mem.Allocator) !void {
    const first = try allocator.alloc(u8, 16);
    const second = allocator.alloc(u8, 32) catch |err| return err;
    allocator.free(second);
    allocator.free(first);
}

fn swallowsTheFailure(allocator: std.mem.Allocator) !void {
    const first = allocator.alloc(u8, 16) catch return;
    allocator.free(first);
    const second = allocator.alloc(u8, 32) catch return;
    allocator.free(second);
}

fn failsWithSomethingElse(allocator: std.mem.Allocator) !void {
    const first = try allocator.alloc(u8, 16);
    defer allocator.free(first);
    return error.Unrelated;
}

/// Allocates three times on its first call and two on every call after, which
/// is the shape of the real defect: an allocation count that varies between
/// runs of the same input. Nothing here leaks and nothing is hidden.
fn allocatesFewerTimesAfterTheFirstCall(allocator: std.mem.Allocator, calls: *usize) !void {
    const wanted: usize = if (calls.* == 0) 3 else 2;
    calls.* += 1;
    var made: [3][]u8 = undefined;
    var count: usize = 0;
    errdefer for (made[0..count]) |block| allocator.free(block);
    while (count < wanted) : (count += 1) {
        made[count] = try allocator.alloc(u8, 8);
    }
    for (made[0..count]) |block| allocator.free(block);
}

test "a well-behaved function passes" {
    try checkAllocationFailures(std.testing.allocator, allocatesTwiceAndCleansUp, .{});
}

test "a leak on the out-of-memory path is still caught" {
    // Backed by an arena rather than the testing allocator: the planted leak is
    // the point of the test, and the testing allocator would report it a second
    // time as the test's own leak. The arena reclaims it at scope exit.
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    try std.testing.expectError(
        error.MemoryLeakDetected,
        checkAllocationFailuresInternal(arena.allocator(), leaksTheFirstAllocation, .{}, false),
    );
}

test "a swallowed out-of-memory error is still caught" {
    try std.testing.expectError(
        error.SwallowedOutOfMemoryError,
        checkAllocationFailures(std.testing.allocator, swallowsTheFailure, .{}),
    );
}

test "an unrelated error propagates rather than being read as a memory finding" {
    try std.testing.expectError(
        error.Unrelated,
        checkAllocationFailures(std.testing.allocator, failsWithSomethingElse, .{}),
    );
}

test "an allocation count that varies after in-place growth is disabled is rejected" {
    var calls: usize = 0;
    try std.testing.expectError(
        error.NondeterministicMemoryUsage,
        checkAllocationFailures(
            std.testing.allocator,
            allocatesFewerTimesAfterTheFirstCall,
            .{&calls},
        ),
    );
}
