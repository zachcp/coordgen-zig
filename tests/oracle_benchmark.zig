//! Absolute oracle timing baseline for the versioned representative corpus.
//!
//! This intentionally does not claim native-vs-oracle parity: there is no
//! native generation entry point yet. It establishes the executable half of
//! the frozen method so the eventual native implementation is timed in this
//! same process, build, target, optimize mode, and member loop.

const std = @import("std");
const builtin = @import("builtin");
const conformance = @import("conformance");
const c_abi = @import("c_abi");

const corpus = conformance.corpus;
const repetitions = 3;

extern fn coordgen_generate(input: *const c_abi.Input, result: *c_abi.Result) u32;
extern fn coordgen_result_free(result: *c_abi.Result) void;
extern fn cgz_benchmark_now_ns() u64;

const Member = struct {
    partition: corpus.Partition,
    index: u32,
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
    var samples: [5]std.ArrayList(u64) = @splat(.empty);
    defer for (&samples) |*bucket| bucket.deinit(gpa);

    for (members) |selected| {
        const molecule = try corpus.generate(gpa, selected.partition, selected.index);
        defer molecule.deinit(gpa);
        const bucket = corpus.SizeBucket.of(molecule.atoms.len);
        const prepared = try PreparedInput.init(gpa, molecule);
        defer prepared.deinit(gpa);

        // One untimed generation initializes any process-local cold state.
        try generateOnce(&prepared.input);
        for (0..repetitions) |_| {
            const start = cgz_benchmark_now_ns();
            if (start == 0) return error.MonotonicClockUnavailable;
            try generateOnce(&prepared.input);
            const finish = cgz_benchmark_now_ns();
            if (finish <= start) return error.MonotonicClockUnavailable;
            try samples[@intFromEnum(bucket)].append(gpa, finish - start);
        }
    }

    var file = std.Io.File.stdout();
    var buffer: [4096]u8 = undefined;
    var file_writer = file.writer(io, &buffer);
    const out = &file_writer.interface;
    try out.print("# target\t{t}-{t}\n", .{ builtin.target.cpu.arch, builtin.target.os.tag });
    try out.print("# toolchain\t{s}\n", .{builtin.zig_version_string});
    try out.print("# optimize\t{t}\n", .{builtin.mode});
    try out.writeAll("# implementation\tbucket\tmembers\tsamples\tmedian_ns\tp95_ns\n");
    for (&samples, 0..) |*bucket_samples, raw_bucket| {
        if (bucket_samples.items.len != 20 * repetitions) {
            return error.InvalidRepresentativePopulation;
        }
        std.mem.sort(u64, bucket_samples.items, {}, u64LessThan);
        const median = percentile(bucket_samples.items, 50);
        const p95 = percentile(bucket_samples.items, 95);
        const bucket: corpus.SizeBucket = @enumFromInt(raw_bucket);
        try out.print("oracle\t{t}\t20\t{d}\t{d}\t{d}\n", .{
            bucket,
            bucket_samples.items.len,
            median,
            p95,
        });
    }
    try out.flush();
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

fn generateOnce(input: *const c_abi.Input) !void {
    var result: c_abi.Result = .{};
    const status = coordgen_generate(input, &result);
    if (status != 0) return error.OracleGenerationFailed;
    coordgen_result_free(&result);
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

test "nearest-rank percentiles do not interpolate timings" {
    const values = [_]u64{ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
    try std.testing.expectEqual(@as(u64, 5), percentile(&values, 50));
    try std.testing.expectEqual(@as(u64, 10), percentile(&values, 95));
}
