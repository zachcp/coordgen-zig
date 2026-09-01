//! Native-vs-oracle timing on the versioned representative corpus
//! (cgz-7v2.4.7).
//!
//! Both implementations are timed in ONE process, on the same members, in the
//! same build, target and optimize mode, inside the same member loop - the
//! identity rule `RunIdentity.requireSameBuild()` already enforces for
//! coordinates, applied to timings.
//!
//! ## One asymmetry, stated rather than hidden
//!
//! The oracle is reached through the C entry points it provides; native is
//! reached through its Zig API. That is not a preference: the linked oracle
//! already defines the `coordgen_*` symbols, and a second definition is the
//! cgz-r28 duplicate-symbol defect, so native's own C layer cannot be in this
//! binary. The difference is native's DTO conversion - a bounded copy of the
//! atom and bond arrays, with no layout work in it - which is charged to the
//! ORACLE side here and not to native. The ratios below are therefore mildly
//! generous to native, and a threshold set from them inherits that. Anyone
//! tightening these numbers should close that gap first.
//!
//! With `--enforce`, the same binary compares each bucket's ratio against
//! conformance/performance_thresholds.tsv and exits non-zero on a regression,
//! so the gate reads exactly the numbers the baseline printed.

const std = @import("std");
const builtin = @import("builtin");
const conformance = @import("conformance");
const c_abi = @import("c_abi");
const api = @import("api");

const thresholds_table = @embedFile("performance_thresholds");

const corpus = conformance.corpus;
const repetitions = 3;

extern fn coordgen_generate(input: *const c_abi.Input, result: *c_abi.Result) u32;
extern fn coordgen_result_free(result: *c_abi.Result) void;
extern fn cgz_benchmark_now_ns() u64;
extern fn cgz_benchmark_thread_cpu_ns() u64;

const Member = struct {
    partition: corpus.Partition,
    index: u32,
};

const MemberTiming = struct {
    member: Member,
    bucket: corpus.SizeBucket,
    oracle_ns: [repetitions]u64,
    oracle_cpu_ns: [repetitions]u64,
    native_ns: [repetitions]u64,
    native_cpu_ns: [repetitions]u64,
};

// Twenty members per atom-count bucket. The drug-like population is included
// in full; remaining slots are the lowest-index adversarial members in that
// bucket. Changing this reviewed population changes the benchmark and requires
// a decision Bead rather than silently making a regression disappear.
const members = [_]Member{
    .{ .partition = .drug_like, .index = 5 },
    .{ .partition = .adversarial, .index = 1 }, .{ .partition = .adversarial, .index = 8 },
    .{ .partition = .adversarial, .index = 11 }, .{ .partition = .adversarial, .index = 30 },
    .{ .partition = .adversarial, .index = 32 }, .{ .partition = .adversarial, .index = 42 },
    .{ .partition = .adversarial, .index = 59 }, .{ .partition = .adversarial, .index = 61 },
    .{ .partition = .adversarial, .index = 70 }, .{ .partition = .adversarial, .index = 78 },
    .{ .partition = .adversarial, .index = 81 }, .{ .partition = .adversarial, .index = 94 },
    .{ .partition = .adversarial, .index = 108 }, .{ .partition = .adversarial, .index = 109 },
    .{ .partition = .adversarial, .index = 111 }, .{ .partition = .adversarial, .index = 118 },
    .{ .partition = .adversarial, .index = 138 }, .{ .partition = .adversarial, .index = 168 },
    .{ .partition = .adversarial, .index = 198 },

    .{ .partition = .drug_like, .index = 0 }, .{ .partition = .drug_like, .index = 1 },
    .{ .partition = .drug_like, .index = 2 }, .{ .partition = .drug_like, .index = 3 },
    .{ .partition = .drug_like, .index = 4 }, .{ .partition = .drug_like, .index = 6 },
    .{ .partition = .adversarial, .index = 4 }, .{ .partition = .adversarial, .index = 5 },
    .{ .partition = .adversarial, .index = 7 }, .{ .partition = .adversarial, .index = 9 },
    .{ .partition = .adversarial, .index = 13 }, .{ .partition = .adversarial, .index = 17 },
    .{ .partition = .adversarial, .index = 19 }, .{ .partition = .adversarial, .index = 20 },
    .{ .partition = .adversarial, .index = 22 }, .{ .partition = .adversarial, .index = 23 },
    .{ .partition = .adversarial, .index = 28 }, .{ .partition = .adversarial, .index = 33 },
    .{ .partition = .adversarial, .index = 35 }, .{ .partition = .adversarial, .index = 36 },

    .{ .partition = .adversarial, .index = 2 }, .{ .partition = .adversarial, .index = 6 },
    .{ .partition = .adversarial, .index = 10 }, .{ .partition = .adversarial, .index = 12 },
    .{ .partition = .adversarial, .index = 25 }, .{ .partition = .adversarial, .index = 26 },
    .{ .partition = .adversarial, .index = 27 }, .{ .partition = .adversarial, .index = 29 },
    .{ .partition = .adversarial, .index = 34 }, .{ .partition = .adversarial, .index = 41 },
    .{ .partition = .adversarial, .index = 46 }, .{ .partition = .adversarial, .index = 49 },
    .{ .partition = .adversarial, .index = 53 }, .{ .partition = .adversarial, .index = 55 },
    .{ .partition = .adversarial, .index = 63 }, .{ .partition = .adversarial, .index = 64 },
    .{ .partition = .adversarial, .index = 65 }, .{ .partition = .adversarial, .index = 66 },
    .{ .partition = .adversarial, .index = 67 }, .{ .partition = .adversarial, .index = 68 },

    .{ .partition = .adversarial, .index = 0 }, .{ .partition = .adversarial, .index = 3 },
    .{ .partition = .adversarial, .index = 14 }, .{ .partition = .adversarial, .index = 15 },
    .{ .partition = .adversarial, .index = 16 }, .{ .partition = .adversarial, .index = 18 },
    .{ .partition = .adversarial, .index = 24 }, .{ .partition = .adversarial, .index = 31 },
    .{ .partition = .adversarial, .index = 37 }, .{ .partition = .adversarial, .index = 38 },
    .{ .partition = .adversarial, .index = 43 }, .{ .partition = .adversarial, .index = 45 },
    .{ .partition = .adversarial, .index = 51 }, .{ .partition = .adversarial, .index = 56 },
    .{ .partition = .adversarial, .index = 57 }, .{ .partition = .adversarial, .index = 62 },
    .{ .partition = .adversarial, .index = 71 }, .{ .partition = .adversarial, .index = 73 },
    .{ .partition = .adversarial, .index = 77 }, .{ .partition = .adversarial, .index = 82 },

    .{ .partition = .adversarial, .index = 21 }, .{ .partition = .adversarial, .index = 139 },
    .{ .partition = .adversarial, .index = 154 }, .{ .partition = .adversarial, .index = 188 },
    .{ .partition = .adversarial, .index = 242 }, .{ .partition = .adversarial, .index = 264 },
    .{ .partition = .adversarial, .index = 292 }, .{ .partition = .adversarial, .index = 301 },
    .{ .partition = .adversarial, .index = 340 }, .{ .partition = .adversarial, .index = 360 },
    .{ .partition = .adversarial, .index = 380 }, .{ .partition = .adversarial, .index = 458 },
    .{ .partition = .adversarial, .index = 482 }, .{ .partition = .adversarial, .index = 491 },
    .{ .partition = .adversarial, .index = 540 }, .{ .partition = .adversarial, .index = 570 },
    .{ .partition = .adversarial, .index = 596 }, .{ .partition = .adversarial, .index = 752 },
    .{ .partition = .adversarial, .index = 780 }, .{ .partition = .adversarial, .index = 791 },
};

pub fn main(init: std.process.Init) !void {
    const gpa = init.gpa;
    const io = init.io;
    var enforce = false;
    var arguments = init.minimal.args.iterate();
    _ = arguments.next();
    while (arguments.next()) |argument| {
        if (std.mem.eql(u8, argument, "--enforce")) {
            enforce = true;
        } else {
            return error.UnknownArgument;
        }
    }

    var samples: [5]std.ArrayList(u64) = @splat(.empty);
    defer for (&samples) |*bucket| bucket.deinit(gpa);
    var native_samples: [5]std.ArrayList(u64) = @splat(.empty);
    defer for (&native_samples) |*bucket| bucket.deinit(gpa);
    var member_timings: std.ArrayList(MemberTiming) = .empty;
    defer member_timings.deinit(gpa);
    // Members native refused, per bucket. Reported, because a bucket where
    // native handles three of twenty members has a ratio that says very
    // little, and the reader has to be able to see that.
    var unsupported: [5]usize = @splat(0);

    for (members, 0..) |selected, member_ordinal| {
        const molecule = try corpus.generate(gpa, selected.partition, selected.index);
        defer molecule.deinit(gpa);
        const bucket = corpus.SizeBucket.of(molecule.atoms.len);
        const prepared = try PreparedInput.init(gpa, molecule);
        defer prepared.deinit(gpa);
        const native = try NativeInput.init(gpa, molecule);
        defer native.deinit(gpa);

        // Native still refuses domains the oracle accepts. A member it cannot
        // lay out is SKIPPED ON BOTH SIDES and counted, rather than timed for
        // the oracle only: a ratio between medians taken over different member
        // subsets is not a like-for-like ratio, and dropping native's hard
        // members while keeping the oracle's would flatter native by exactly
        // the amount that matters.
        const native_supported = nativeGenerateOnce(gpa, native.input()) catch |err| switch (err) {
            error.EmptyGraph,
            error.TooManyItems,
            error.InvalidAtomicNumber,
            error.InvalidBondOrder,
            error.InvalidAtomIndex,
            error.InvalidStereo,
            error.InvalidCoordinate,
            error.InvalidOption,
            error.InvalidMapping,
            error.Unsupported,
            => false,
            error.OutOfMemory => return err,
        };
        if (!native_supported) {
            unsupported[@intFromEnum(bucket)] += 1;
            continue;
        }
        // One untimed generation of each initializes process-local cold state.
        try generateOnce(&prepared.input);
        var timing: MemberTiming = .{
            .member = selected,
            .bucket = bucket,
            .oracle_ns = undefined,
            .oracle_cpu_ns = undefined,
            .native_ns = undefined,
            .native_cpu_ns = undefined,
        };
        // Alternate which implementation runs first. The previous grouped
        // oracle-then-native order made the oracle an invalid control for a
        // quota, scheduling, or thermal effect that began during the native
        // half of a member. Vary the first implementation by both repetition
        // and member so neither side systematically consumes a fresh interval.
        for (0..repetitions) |repetition| {
            if (oracleRunsFirst(member_ordinal, repetition)) {
                const oracle_elapsed = try timeOracle(&prepared.input);
                timing.oracle_ns[repetition] = oracle_elapsed.wall_ns;
                timing.oracle_cpu_ns[repetition] = oracle_elapsed.cpu_ns;
                const native_elapsed = try timeNative(gpa, native.input());
                timing.native_ns[repetition] = native_elapsed.wall_ns;
                timing.native_cpu_ns[repetition] = native_elapsed.cpu_ns;
            } else {
                const native_elapsed = try timeNative(gpa, native.input());
                timing.native_ns[repetition] = native_elapsed.wall_ns;
                timing.native_cpu_ns[repetition] = native_elapsed.cpu_ns;
                const oracle_elapsed = try timeOracle(&prepared.input);
                timing.oracle_ns[repetition] = oracle_elapsed.wall_ns;
                timing.oracle_cpu_ns[repetition] = oracle_elapsed.cpu_ns;
            }
            try samples[@intFromEnum(bucket)].append(gpa, timing.oracle_ns[repetition]);
            try native_samples[@intFromEnum(bucket)].append(gpa, timing.native_ns[repetition]);
        }
        try member_timings.append(gpa, timing);
    }

    var file = std.Io.File.stdout();
    var buffer: [4096]u8 = undefined;
    var file_writer = file.writer(io, &buffer);
    const out = &file_writer.interface;
    try out.print("# target\t{t}-{t}\n", .{ builtin.target.cpu.arch, builtin.target.os.tag });
    try out.print("# toolchain\t{s}\n", .{builtin.zig_version_string});
    try out.print("# optimize\t{t}\n", .{builtin.mode});
    try out.writeAll(
        "# ratio rows carry compared/unsupported in the members/samples columns;" ++
            " a member native cannot lay out is skipped on BOTH sides so the" ++
            " ratio stays like-for-like\n",
    );
    try out.writeAll(
        "# member timings are emitted after measurement so diagnostics do not perturb later samples\n" ++
            "# member\tpartition\tindex\tbucket",
    );
    for (0..repetitions) |repetition| try out.print("\toracle_ns_{d}", .{repetition});
    for (0..repetitions) |repetition| try out.print("\toracle_cpu_ns_{d}", .{repetition});
    for (0..repetitions) |repetition| try out.print("\tnative_ns_{d}", .{repetition});
    for (0..repetitions) |repetition| try out.print("\tnative_cpu_ns_{d}", .{repetition});
    try out.writeByte('\n');
    for (member_timings.items) |timing| {
        try out.print("# member\t{t}\t{d}\t{t}", .{
            timing.member.partition,
            timing.member.index,
            timing.bucket,
        });
        for (timing.oracle_ns) |sample| try out.print("\t{d}", .{sample});
        for (timing.oracle_cpu_ns) |sample| try out.print("\t{d}", .{sample});
        for (timing.native_ns) |sample| try out.print("\t{d}", .{sample});
        for (timing.native_cpu_ns) |sample| try out.print("\t{d}", .{sample});
        try out.writeByte('\n');
    }
    try out.writeAll("# implementation\tbucket\tmembers\tsamples\tmedian_ns\tp95_ns\n");
    var violations: usize = 0;
    for (&samples, &native_samples, 0..) |*bucket_samples, *bucket_native, raw_bucket| {
        const compared = 20 - unsupported[raw_bucket];
        if (bucket_samples.items.len != compared * repetitions) {
            return error.InvalidRepresentativePopulation;
        }
        if (bucket_samples.items.len == 0) {
            const empty_bucket: corpus.SizeBucket = @enumFromInt(raw_bucket);
            try out.print("oracle\t{t}\t0\t0\t-\t-\t(native supports none of 20)\n", .{empty_bucket});
            if (enforce) {
                // A bucket that used to be comparable and is not any more is a
                // domain regression, and it would otherwise be invisible: with
                // no members timed there is no ratio to exceed.
                const limits = thresholdFor(empty_bucket) orelse {
                    try out.print(
                        "FAIL\t{t}\tno row in conformance/performance_thresholds.tsv\n",
                        .{empty_bucket},
                    );
                    violations += 1;
                    continue;
                };
                if (limits.min_compared != 0) {
                    try out.print(
                        "FAIL\t{t}\tnative laid out 0 members; the recorded floor is {d}\n",
                        .{ empty_bucket, limits.min_compared },
                    );
                    violations += 1;
                }
            }
            continue;
        }
        std.mem.sort(u64, bucket_samples.items, {}, u64LessThan);
        const median = percentile(bucket_samples.items, 50);
        const p95 = percentile(bucket_samples.items, 95);
        const bucket: corpus.SizeBucket = @enumFromInt(raw_bucket);
        try out.print("oracle\t{t}\t{d}\t{d}\t{d}\t{d}\n", .{
            bucket,
            compared,
            bucket_samples.items.len,
            median,
            p95,
        });
        std.mem.sort(u64, bucket_native.items, {}, u64LessThan);
        const native_median = percentile(bucket_native.items, 50);
        const native_p95 = percentile(bucket_native.items, 95);
        try out.print("native\t{t}\t{d}\t{d}\t{d}\t{d}\n", .{
            bucket,
            bucket_native.items.len / repetitions,
            bucket_native.items.len,
            native_median,
            native_p95,
        });
        const median_ratio = ratioOf(native_median, median);
        const p95_ratio = ratioOf(native_p95, p95);
        try out.print("ratio\t{t}\t{d}\t{d}\t{d:.3}\t{d:.3}\n", .{
            bucket,
            compared,
            unsupported[raw_bucket],
            median_ratio,
            p95_ratio,
        });
        if (!enforce) continue;
        const limits = thresholdFor(bucket) orelse {
            try out.print(
                "FAIL\t{t}\tno row in conformance/performance_thresholds.tsv\n",
                .{bucket},
            );
            violations += 1;
            continue;
        };
        if (compared < limits.min_compared) {
            try out.print(
                "FAIL\t{t}\tnative laid out {d} of 20 members; the recorded floor is {d}\n",
                .{ bucket, compared, limits.min_compared },
            );
            violations += 1;
        }
        if (median_ratio > limits.median) {
            try out.print("FAIL\t{t}\tmedian ratio {d:.3} exceeds {d:.3}\n", .{
                bucket, median_ratio, limits.median,
            });
            violations += 1;
        }
        if (p95_ratio > limits.p95) {
            try out.print("FAIL\t{t}\tp95 ratio {d:.3} exceeds {d:.3}\n", .{
                bucket, p95_ratio, limits.p95,
            });
            violations += 1;
        }
    }
    try out.flush();
    if (enforce and violations != 0) return error.PerformanceRegression;
}

fn ratioOf(native: u64, oracle: u64) f64 {
    if (oracle == 0) return std.math.inf(f64);
    return @as(f64, @floatFromInt(native)) / @as(f64, @floatFromInt(oracle));
}

fn oracleRunsFirst(member_ordinal: usize, repetition: usize) bool {
    return (member_ordinal + repetition) % 2 == 0;
}

const Threshold = struct { median: f64, p95: f64, min_compared: usize };

/// The committed per-bucket limits. Absent row means absent threshold, which
/// `--enforce` treats as a failure rather than as permission.
fn thresholdFor(bucket: corpus.SizeBucket) ?Threshold {
    var lines = std.mem.splitScalar(u8, thresholds_table, '\n');
    while (lines.next()) |line| {
        if (line.len == 0 or line[0] == '#') continue;
        var fields = std.mem.splitScalar(u8, line, '\t');
        const name = fields.next() orelse continue;
        const median = fields.next() orelse continue;
        const p95 = fields.next() orelse continue;
        if (!std.mem.eql(u8, name, @tagName(bucket))) continue;
        const compared = fields.next() orelse return null;
        return .{
            .median = std.fmt.parseFloat(f64, std.mem.trim(u8, median, " \r")) catch return null,
            .p95 = std.fmt.parseFloat(f64, std.mem.trim(u8, p95, " \r")) catch return null,
            .min_compared = std.fmt.parseInt(usize, std.mem.trim(u8, compared, " \r"), 10) catch return null,
        };
    }
    return null;
}

const PreparedInput = struct {
    atoms: []c_abi.AtomInput,
    bonds: []c_abi.BondInput,
    input: c_abi.Input,

    fn init(gpa: std.mem.Allocator, molecule: corpus.Molecule) !PreparedInput {
        const atoms = try gpa.alloc(c_abi.AtomInput, molecule.atoms.len);
        errdefer gpa.free(atoms);
        for (atoms, molecule.atoms) |*atom, source| {
            atom.* = .{ .atomic_number = source.atomic_number, .formal_charge = source.formal_charge };
        }
        const bonds = try gpa.alloc(c_abi.BondInput, molecule.bonds.len);
        errdefer gpa.free(bonds);
        for (bonds, molecule.bonds) |*bond, source| {
            bond.* = .{ .start = source.start, .end = source.end, .order = source.order };
        }
        return .{
            .atoms = atoms,
            .bonds = bonds,
            .input = .{
                .atoms = .{ .ptr = atoms.ptr, .len = @intCast(atoms.len) },
                .bonds = .{ .ptr = bonds.ptr, .len = @intCast(bonds.len) },
            },
        };
    }

    fn deinit(self: PreparedInput, gpa: std.mem.Allocator) void {
        gpa.free(self.bonds);
        gpa.free(self.atoms);
    }
};

/// The same members, shaped for the Zig API. Built once per member and reused
/// across repetitions, exactly as `PreparedInput` is, so neither side is
/// charged for construction.
const NativeInput = struct {
    atoms: []api.AtomInput,
    bonds: []api.BondInput,

    fn init(gpa: std.mem.Allocator, molecule: corpus.Molecule) !NativeInput {
        const atoms = try gpa.alloc(api.AtomInput, molecule.atoms.len);
        errdefer gpa.free(atoms);
        for (atoms, molecule.atoms) |*atom, source| {
            atom.* = .{
                .atomic_number = @enumFromInt(source.atomic_number),
                .formal_charge = source.formal_charge,
            };
        }
        const bonds = try gpa.alloc(api.BondInput, molecule.bonds.len);
        errdefer gpa.free(bonds);
        for (bonds, molecule.bonds) |*bond, source| {
            bond.* = .{
                .start = source.start,
                .end = source.end,
                .order = @enumFromInt(source.order),
            };
        }
        return .{ .atoms = atoms, .bonds = bonds };
    }

    fn input(self: NativeInput) api.Input {
        return .{ .atoms = self.atoms, .bonds = self.bonds };
    }

    fn deinit(self: NativeInput, gpa: std.mem.Allocator) void {
        gpa.free(self.bonds);
        gpa.free(self.atoms);
    }
};

fn nativeGenerateOnce(gpa: std.mem.Allocator, input: api.Input) !bool {
    var result = try api.generate(gpa, input);
    result.deinit();
    return true;
}

const Elapsed = struct { wall_ns: u64, cpu_ns: u64 };

fn timeNative(gpa: std.mem.Allocator, input: api.Input) !Elapsed {
    const wall_start = cgz_benchmark_now_ns();
    const cpu_start = cgz_benchmark_thread_cpu_ns();
    if (wall_start == 0 or cpu_start == 0) return error.MonotonicClockUnavailable;
    _ = try nativeGenerateOnce(gpa, input);
    const cpu_finish = cgz_benchmark_thread_cpu_ns();
    const wall_finish = cgz_benchmark_now_ns();
    if (wall_finish <= wall_start or cpu_finish <= cpu_start) return error.MonotonicClockUnavailable;
    return .{ .wall_ns = wall_finish - wall_start, .cpu_ns = cpu_finish - cpu_start };
}

fn generateOnce(input: *const c_abi.Input) !void {
    var result: c_abi.Result = .{};
    const status = coordgen_generate(input, &result);
    if (status != 0) return error.OracleGenerationFailed;
    coordgen_result_free(&result);
}

fn timeOracle(input: *const c_abi.Input) !Elapsed {
    const wall_start = cgz_benchmark_now_ns();
    const cpu_start = cgz_benchmark_thread_cpu_ns();
    if (wall_start == 0 or cpu_start == 0) return error.MonotonicClockUnavailable;
    try generateOnce(input);
    const cpu_finish = cgz_benchmark_thread_cpu_ns();
    const wall_finish = cgz_benchmark_now_ns();
    if (wall_finish <= wall_start or cpu_finish <= cpu_start) return error.MonotonicClockUnavailable;
    return .{ .wall_ns = wall_finish - wall_start, .cpu_ns = cpu_finish - cpu_start };
}

fn percentile(sorted: []const u64, percent: usize) u64 {
    std.debug.assert(sorted.len > 0 and percent > 0 and percent <= 100);
    return sorted[((sorted.len * percent + 99) / 100) - 1];
}

fn u64LessThan(_: void, left: u64, right: u64) bool {
    return left < right;
}

test "representative benchmark has twenty members in every size bucket" {
    var counts: [5]usize = @splat(0);
    for (members) |selected| {
        const molecule = try corpus.generate(std.testing.allocator, selected.partition, selected.index);
        defer molecule.deinit(std.testing.allocator);
        counts[@intFromEnum(corpus.SizeBucket.of(molecule.atoms.len))] += 1;
    }
    try std.testing.expectEqual([_]usize{ 20, 20, 20, 20, 20 }, counts);
}

test "timing order is balanced and changes within every member" {
    var oracle_first: usize = 0;
    for (members, 0..) |_, member_ordinal| {
        for (0..repetitions) |repetition| {
            if (oracleRunsFirst(member_ordinal, repetition)) oracle_first += 1;
            if (repetition != 0) {
                try std.testing.expect(oracleRunsFirst(member_ordinal, repetition) !=
                    oracleRunsFirst(member_ordinal, repetition - 1));
            }
        }
    }
    try std.testing.expectEqual(members.len * repetitions / 2, oracle_first);
}

test "nearest-rank percentiles do not interpolate timings" {
    const values = [_]u64{ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
    try std.testing.expectEqual(@as(u64, 5), percentile(&values, 50));
    try std.testing.expectEqual(@as(u64, 10), percentile(&values, 95));
}
