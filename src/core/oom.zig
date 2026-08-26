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

/// Count the DISTINCT allocation sites a call reaches, where a site is the
/// return address the allocator is called through.
///
/// This is a different quantity from the number of allocation indices
/// `checkAllocationFailures` sweeps, and it answers a different question. The
/// sweep proves every index it reaches is handled; this proves the SET of
/// sites reached does not silently shrink as the code grows. A representative
/// input that stops reaching a module still passes every sweep - it just
/// stops testing that module - and nothing else in the tree would notice.
/// docs/architecture/SUCCESS_CRITERIA.md requires it for exactly that reason.
///
/// Counted under the same no-in-place allocator the sweep uses, so a growable
/// container's fallback path is reached deterministically rather than
/// depending on the residual heap.
///
/// The count is a property of (source, toolchain, optimize mode), because
/// inlining decides how many machine call sites one source call site becomes.
/// It is comparable against a recorded count taken under the same three, which
/// is how the committed table is keyed. The addresses themselves are never
/// recorded: they are not stable across builds and would be meaningless in a
/// diff.
pub fn countAllocationSites(
    backing_allocator: std.mem.Allocator,
    comptime test_fn: anytype,
    extra_args: anytype,
) !usize {
    var no_in_place = NoInPlaceAllocator{ .backing = backing_allocator };
    var counter = SiteCounter{ .backing = no_in_place.allocator(), .sites = .empty };
    defer counter.sites.deinit(backing_allocator);
    counter.owner = backing_allocator;
    try @call(.auto, test_fn, .{counter.allocator()} ++ extra_args);
    if (counter.overflow) return error.OutOfMemory;
    return counter.sites.count();
}

/// Records the return address of every allocation, and otherwise forwards.
/// Only `alloc` is recorded: `resize` and `remap` are declined by the wrapper
/// underneath, so a growth reaches this as a fresh `alloc` at the same site.
const SiteCounter = struct {
    backing: std.mem.Allocator,
    owner: std.mem.Allocator = undefined,
    sites: std.AutoHashMapUnmanaged(usize, void),
    /// Set when the site set itself could not grow. Reported as an error
    /// rather than silently undercounting.
    overflow: bool = false,

    fn allocator(self: *SiteCounter) std.mem.Allocator {
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
        const self: *SiteCounter = @ptrCast(@alignCast(context));
        const result = self.backing.rawAlloc(len, alignment, return_address);
        if (result != null) {
            self.sites.put(self.owner, return_address, {}) catch {
                self.overflow = true;
            };
        }
        return result;
    }

    fn resize(context: *anyopaque, memory: []u8, alignment: std.mem.Alignment, new_len: usize, return_address: usize) bool {
        const self: *SiteCounter = @ptrCast(@alignCast(context));
        return self.backing.rawResize(memory, alignment, new_len, return_address);
    }

    fn remap(context: *anyopaque, memory: []u8, alignment: std.mem.Alignment, new_len: usize, return_address: usize) ?[*]u8 {
        const self: *SiteCounter = @ptrCast(@alignCast(context));
        return self.backing.rawRemap(memory, alignment, new_len, return_address);
    }

    fn free(context: *anyopaque, memory: []u8, alignment: std.mem.Alignment, return_address: usize) void {
        const self: *SiteCounter = @ptrCast(@alignCast(context));
        self.backing.rawFree(memory, alignment, return_address);
    }
};

/// Look up one row of the committed allocation-site table. The table's bytes
/// are passed in rather than embedded here: `@embedFile` resolves against the
/// importing module's root, and the two suites that need this table live in
/// different modules. Keeping the PARSER single-sourced is the part that
/// matters - two copies of it drifting apart is the duplication cgz-7v2.24
/// spent its time removing.
///
/// Returns null when no row matches, which callers must treat as a failure
/// rather than a skip: a mode with no committed rows would otherwise pass
/// vacuously (cgz-r20).
pub fn recordedSiteFloor(
    table: []const u8,
    entry_point: []const u8,
    member: []const u8,
    optimize_mode: []const u8,
) ?usize {
    var lines = std.mem.splitScalar(u8, table, '\n');
    while (lines.next()) |line| {
        if (line.len == 0 or line[0] == '#') continue;
        var fields = std.mem.splitScalar(u8, line, '\t');
        const row_entry = fields.next() orelse continue;
        const row_member = fields.next() orelse continue;
        const row_mode = fields.next() orelse continue;
        const row_value = fields.next() orelse continue;
        if (!std.mem.eql(u8, row_entry, entry_point)) continue;
        if (!std.mem.eql(u8, row_member, member)) continue;
        if (!std.mem.eql(u8, row_mode, optimize_mode)) continue;
        return std.fmt.parseInt(usize, std.mem.trim(u8, row_value, " \r"), 10) catch null;
    }
    return null;
}

/// The name this build's optimize mode is recorded under. The count is not
/// comparable across modes - inlining decides how many machine call sites one
/// source call site becomes - so every row is keyed by it.
pub const optimize_mode_name: []const u8 = switch (@import("builtin").mode) {
    .Debug => "debug",
    .ReleaseSafe => "release_safe",
    .ReleaseFast => "release_fast",
    .ReleaseSmall => "release_small",
};

/// Assert that `measured` has not fallen below the table's floor, reporting by
/// name with enough context to act on.
pub fn expectSiteFloor(
    table: []const u8,
    entry_point: []const u8,
    member: []const u8,
    measured: usize,
) !void {
    const mode = optimize_mode_name;
    const floor = recordedSiteFloor(table, entry_point, member, mode) orelse {
        std.debug.print(
            "conformance/allocation_sites.tsv has no row for {s}/{s}/{s};" ++
                " measure and commit it rather than skipping the gate\n",
            .{ entry_point, member, mode },
        );
        return error.NoRecordedAllocationSites;
    };
    if (measured < floor) {
        std.debug.print(
            "allocation site coverage dropped for {s}/{s}/{s}: {d} sites," ++
                " recorded floor {d}. An input that stops reaching a module still" ++
                " passes every allocation sweep - it just stops testing that module.\n",
            .{ entry_point, member, mode, measured, floor },
        );
        return error.AllocationSiteCoverageDecreased;
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

fn allocatesAtTwoSitesRepeatedly(allocator: std.mem.Allocator) !void {
    var first: [3][]u8 = undefined;
    for (&first) |*block| block.* = try allocator.alloc(u8, 16);
    defer for (first) |block| allocator.free(block);
    var second: [2][]u8 = undefined;
    for (&second) |*block| block.* = try allocator.alloc(u8, 32);
    defer for (second) |block| allocator.free(block);
}

test "site counting counts sites, not allocations" {
    // Five allocations from two source sites. A count of five would mean the
    // helper is measuring the same thing checkAllocationFailures already does.
    //
    // The exact number is deliberately not asserted: it is 2 in Debug and 1 in
    // ReleaseSafe, because the two source sites are in one loop body each and
    // inlining collapses them. That mode dependence is real, and it is why the
    // committed table in conformance/allocation_sites.tsv is keyed by optimize
    // mode. What holds in every mode is the separation being claimed here -
    // fewer sites than allocations, and at least one.
    const allocation_count = 5;
    const repeated = try countAllocationSites(std.testing.allocator, allocatesAtTwoSitesRepeatedly, .{});
    try std.testing.expect(repeated >= 1);
    try std.testing.expect(repeated < allocation_count);

    // Two allocations from two sites: nothing to collapse, so this one cannot
    // be below one either way.
    const distinct = try countAllocationSites(std.testing.allocator, allocatesTwiceAndCleansUp, .{});
    try std.testing.expect(distinct >= 1);
    try std.testing.expect(distinct <= 2);
}

test "a missing row fails by name rather than skipping" {
    const table = "# comment\napi.generate\tchain\tdebug\t83\n";
    try std.testing.expectEqual(@as(?usize, 83), recordedSiteFloor(table, "api.generate", "chain", "debug"));
    try std.testing.expectEqual(@as(?usize, null), recordedSiteFloor(table, "api.generate", "chain", "release_fast"));
    try std.testing.expectError(
        error.NoRecordedAllocationSites,
        expectSiteFloor(table, "api.generate", "nonexistent member", 1000),
    );
}

test "a decrease below the recorded floor is caught" {
    const table = "api.generate\tchain\t" ++ optimize_mode_name ++ "\t83\n";
    try expectSiteFloor(table, "api.generate", "chain", 83);
    try expectSiteFloor(table, "api.generate", "chain", 84);
    try std.testing.expectError(
        error.AllocationSiteCoverageDecreased,
        expectSiteFloor(table, "api.generate", "chain", 82),
    );
}
