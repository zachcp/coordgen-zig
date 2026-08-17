//! The native-vs-oracle differential runner (cgz-7v2.4.2).
//!
//! Both implementations run over the same corpus member in one process,
//! build, and target, which is what `COMPARISON_SEMANTICS.md` requires and
//! what makes the architecture axis unobservable here (cgz-r13, cgz-r30).
//!
//! Native is reached through the Zig `api` module, never through the C ABI
//! exports. This binary links the *oracle's* definition of
//! `coordgen_generate`; linking the native library's definition as well is
//! the duplicate-symbol failure cgz-r28 records. The public C entry point has
//! its own end-to-end coverage in `tests/abi_layout.c`.
//!
//! Scope of this runner is the eight public observables. The probe seams
//! (Morgan ranks, rings, fragments, DOFs, template mappings, components) need
//! a native probe surface that does not exist yet; they are named in the
//! parity manifest and remain this bead's outstanding work.
//!
//! Usage: native-oracle-diff --partition NAME [--count N] --ceiling PATH
//!        [--tolerance BOND_LENGTHS] [--baseline PATH]

const std = @import("std");
const api = @import("api");
const conformance = @import("conformance");
const c_abi = @import("c_abi");
const core = @import("core");

const corpus = conformance.corpus;
const comparison = conformance.comparison;

extern fn coordgen_generate(input: *const c_abi.Input, result: *c_abi.Result) u32;
extern fn coordgen_result_free(result: *c_abi.Result) void;

/// The published prior from SUCCESS_CRITERIA.md. It bounds the *oracle's* own
/// float sensitivity, not the port's, and recalibrating it against the first
/// baseline this runner produces is part of cgz-7v2.4.2.
const default_tolerance_bond_lengths: f64 = 0.1;

const public_observables = [_]comparison.Observable{
    .clean_pose,
    .coordinates,
    .input_to_internal,
    .internal_to_input,
    .effective_bond_orders,
    .bond_displays,
    .atom_stereo,
};

const Counter = struct {
    compared: u32 = 0,
    mismatched: u32 = 0,
    /// Pairs the enumeration excused from exact comparison (cgz-r26).
    ceiling_applied: u32 = 0,
};

const Totals = struct {
    members: u32 = 0,
    /// Both implementations returned the same non-ok code. Agreement on a
    /// rejection is agreement, and it is reported rather than hidden.
    agreed_rejections: u32 = 0,
    /// Native declined a domain the oracle produced coordinates for. Never a
    /// silent skip: it is a coverage gap and it fails the run.
    native_unsupported: u32 = 0,
    status_mismatches: u32 = 0,
    max_deviation_bond_lengths: f64 = 0,
    max_deviation_member: [64]u8 = @splat(0),
    max_deviation_member_len: usize = 0,

    fn recordDeviation(self: *Totals, member: []const u8, deviation: f64) void {
        if (deviation <= self.max_deviation_bond_lengths) return;
        self.max_deviation_bond_lengths = deviation;
        const len = @min(member.len, self.max_deviation_member.len);
        @memcpy(self.max_deviation_member[0..len], member[0..len]);
        self.max_deviation_member_len = len;
    }

    fn worstMember(self: *const Totals) []const u8 {
        return self.max_deviation_member[0..self.max_deviation_member_len];
    }
};

pub fn main(init: std.process.Init) !void {
    const arena = init.arena.allocator();
    const gpa = init.gpa;
    const io = init.io;

    var partition: ?corpus.Partition = null;
    var requested_count: u32 = 0;
    var ceiling_path: ?[]const u8 = null;
    var baseline_path: ?[]const u8 = null;
    var tolerance = default_tolerance_bond_lengths;

    const args = try init.minimal.args.toSlice(arena);
    var index: usize = 1;
    while (index < args.len) : (index += 2) {
        if (index + 1 >= args.len) return fatal("option '{s}' needs a value", .{args[index]});
        const value = args[index + 1];
        if (std.mem.eql(u8, args[index], "--partition")) {
            partition = std.meta.stringToEnum(corpus.Partition, value) orelse
                return fatal("unknown partition '{s}'", .{value});
        } else if (std.mem.eql(u8, args[index], "--count")) {
            requested_count = std.fmt.parseInt(u32, value, 10) catch
                return fatal("--count '{s}' is not a number", .{value});
        } else if (std.mem.eql(u8, args[index], "--ceiling")) {
            ceiling_path = value;
        } else if (std.mem.eql(u8, args[index], "--baseline")) {
            baseline_path = value;
        } else if (std.mem.eql(u8, args[index], "--tolerance")) {
            tolerance = std.fmt.parseFloat(f64, value) catch
                return fatal("--tolerance '{s}' is not a number", .{value});
        } else {
            return fatal("unknown option '{s}'", .{args[index]});
        }
    }

    const selected = partition orelse return fatal("--partition is required", .{});
    // Required, not optional: an optional ceiling would let a forgotten
    // argument silently downgrade the run to "no pair is ever excused", which
    // reads as a pass. Same reasoning as corpus-classify's --ceiling.
    const ceiling_file = ceiling_path orelse return fatal("--ceiling is required", .{});
    const count = selected.memberCount(requested_count);
    if (count == 0) return fatal("corpus partition '{t}' would be empty", .{selected});

    const ceiling_text = try std.Io.Dir.cwd().readFileAlloc(io, ceiling_file, arena, .limited(1 << 20));
    const ceiling = comparison.CeilingTable{ .text = ceiling_text };
    ceiling.validateUnique() catch |err|
        return fatal("ceiling file '{s}' is malformed: {t}", .{ ceiling_file, err });

    var report_buffer: [64 * 1024]u8 = undefined;
    var report_file: std.Io.File.Writer = .init(.stderr(), io, &report_buffer);
    const report = &report_file.interface;

    var counters: [public_observables.len]Counter = @splat(.{});
    var totals = Totals{};

    for (0..count) |raw_index| {
        const member_index: u32 = @intCast(raw_index);
        const molecule = try corpus.generate(gpa, selected, member_index);
        defer molecule.deinit(gpa);
        const member = try std.fmt.allocPrint(arena, "{t}/{d}", .{ selected, member_index });
        try compareMember(
            gpa,
            report,
            molecule,
            member,
            ceiling,
            tolerance,
            &counters,
            &totals,
        );
    }

    try writeSummary(report, selected, tolerance, &counters, &totals);
    if (baseline_path) |path| {
        var file = try std.Io.Dir.cwd().createFile(io, path, .{});
        defer file.close(io);
        var buffer: [64 * 1024]u8 = undefined;
        var writer = file.writer(io, &buffer);
        try writeBaseline(&writer.interface, selected, tolerance, &counters, &totals);
        try writer.interface.flush();
    }
    try report.flush();

    if (totals.status_mismatches != 0 or totals.native_unsupported != 0) {
        return error.NativeOracleStatusMismatch;
    }
    for (counters) |counter| {
        if (counter.mismatched != 0) return error.NativeOracleParityMismatch;
    }
}

fn compareMember(
    gpa: std.mem.Allocator,
    report: *std.Io.Writer,
    molecule: corpus.Molecule,
    member: []const u8,
    ceiling: comparison.CeilingTable,
    tolerance: f64,
    counters: *[public_observables.len]Counter,
    totals: *Totals,
) !void {
    totals.members += 1;

    const oracle_atoms = try gpa.alloc(c_abi.AtomInput, molecule.atoms.len);
    defer gpa.free(oracle_atoms);
    const native_atoms = try gpa.alloc(api.AtomInput, molecule.atoms.len);
    defer gpa.free(native_atoms);
    for (oracle_atoms, native_atoms, molecule.atoms) |*oracle_atom, *native_atom, source| {
        oracle_atom.* = .{
            .atomic_number = source.atomic_number,
            .formal_charge = source.formal_charge,
        };
        native_atom.* = .{
            .atomic_number = api.AtomicNumber.fromPublic(source.atomic_number) orelse
                return error.CorpusAtomicNumberOutOfRange,
            .formal_charge = source.formal_charge,
        };
    }

    const oracle_bonds = try gpa.alloc(c_abi.BondInput, molecule.bonds.len);
    defer gpa.free(oracle_bonds);
    const native_bonds = try gpa.alloc(api.BondInput, molecule.bonds.len);
    defer gpa.free(native_bonds);
    for (oracle_bonds, native_bonds, molecule.bonds) |*oracle_bond, *native_bond, source| {
        oracle_bond.* = .{ .start = source.start, .end = source.end, .order = source.order };
        native_bond.* = .{
            .start = source.start,
            .end = source.end,
            .order = api.BondOrder.fromInt(source.order) orelse
                return error.CorpusBondOrderOutOfRange,
        };
    }

    const oracle_input: c_abi.Input = .{
        .atoms = .{ .ptr = oracle_atoms.ptr, .len = @intCast(oracle_atoms.len) },
        .bonds = .{ .ptr = oracle_bonds.ptr, .len = @intCast(oracle_bonds.len) },
    };
    var oracle_result: c_abi.Result = undefined;
    const oracle_status = coordgen_generate(&oracle_input, &oracle_result);
    defer coordgen_result_free(&oracle_result);

    const native_input = api.Input{ .atoms = native_atoms, .bonds = native_bonds };
    var native_value = api.generate(gpa, native_input) catch |err| {
        const native_status = @backingInt(core.errors.code(err));
        if (oracle_status == native_status) {
            totals.agreed_rejections += 1;
            return;
        }
        if (err == error.Unsupported) {
            totals.native_unsupported += 1;
            try report.print(
                "{s}: native declined a domain the oracle accepted (oracle status {d})\n",
                .{ member, oracle_status },
            );
            return;
        }
        totals.status_mismatches += 1;
        try report.print(
            "{s}: status oracle={d} native={d} ({t})\n",
            .{ member, oracle_status, native_status, err },
        );
        return;
    };
    defer native_value.deinit();

    if (oracle_status != 0) {
        totals.status_mismatches += 1;
        try report.print(
            "{s}: status oracle={d} native=0 (native succeeded where the oracle did not)\n",
            .{ member, oracle_status },
        );
        return;
    }

    for (public_observables, 0..) |observable, slot| {
        const row = ceiling.find(member, observable) catch |err|
            return fatal("ceiling lookup failed for {s}: {t}", .{ member, err });
        // A pair carrying a ceiling row is precisely the enumerated
        // order-unstable set (cgz-r26); everything else is held order-stable
        // and therefore to its base tier.
        const plan = comparison.differentialComparison(
            observable,
            row == null,
            row,
            tolerance,
        ) catch |err| {
            try report.print("{s}/{t}: {t}\n", .{ member, observable, err });
            counters[slot].mismatched += 1;
            continue;
        };
        if (row != null) counters[slot].ceiling_applied += 1;
        counters[slot].compared += 1;

        const matched = switch (observable) {
            .coordinates => try compareCoordinateObservable(
                gpa,
                report,
                member,
                oracle_result,
                native_value,
                plan.tolerance_bond_lengths,
                totals,
            ),
            .clean_pose => oracle_result.clean_pose == @intFromBool(native_value.clean_pose),
            .input_to_internal => spansEqual(
                oracle_result.input_to_internal,
                native_value.input_to_internal,
            ),
            .internal_to_input => spansEqual(
                oracle_result.internal_to_input,
                native_value.internal_to_input,
            ),
            .effective_bond_orders => enumSpanEqual(
                oracle_result.effective_bond_orders,
                native_value.effective_bond_orders,
            ),
            .bond_displays => enumSpanEqual(
                oracle_result.bond_displays,
                native_value.bond_displays,
            ),
            .atom_stereo => enumSpanEqual(oracle_result.atom_stereo, native_value.atom_stereo),
            else => unreachable,
        };
        if (!matched) {
            counters[slot].mismatched += 1;
            if (observable != .coordinates) {
                try report.print("{s}/{t}: mismatch\n", .{ member, observable });
            }
        }
    }
}

fn compareCoordinateObservable(
    gpa: std.mem.Allocator,
    report: *std.Io.Writer,
    member: []const u8,
    oracle_result: c_abi.Result,
    native_value: api.Result,
    tolerance: f64,
    totals: *Totals,
) !bool {
    const len: usize = @intCast(oracle_result.coordinates.len);
    if (len != native_value.coordinates.len) {
        try report.print(
            "{s}/coordinates: length oracle={d} native={d}\n",
            .{ member, len, native_value.coordinates.len },
        );
        return false;
    }
    const oracle_points = try gpa.alloc(comparison.Point, len);
    defer gpa.free(oracle_points);
    const native_points = try gpa.alloc(comparison.Point, len);
    defer gpa.free(native_points);
    const raw = oracle_result.coordinates.ptr.?;
    for (oracle_points, native_points, 0..) |*oracle_point, *native_point, i| {
        oracle_point.* = .{ .x = raw[i].x, .y = raw[i].y };
        native_point.* = .{
            .x = native_value.coordinates[i].x,
            .y = native_value.coordinates[i].y,
        };
    }
    // Reflection stays off: a mirrored layout is a different molecule unless
    // the fixture carries an achirality proof (COMPARISON_SEMANTICS.md).
    const result = try comparison.compareCoordinates(
        oracle_points,
        native_points,
        tolerance,
        .{},
    );
    totals.recordDeviation(member, result.max_deviation_bond_lengths);
    if (!result.matches) {
        try report.print(
            "{s}/coordinates: {d:.4} bond lengths exceeds {d:.4}\n",
            .{ member, result.max_deviation_bond_lengths, tolerance },
        );
    }
    return result.matches;
}

fn spansEqual(oracle: c_abi.U32Span, native: []const u32) bool {
    const len: usize = @intCast(oracle.len);
    if (len != native.len) return false;
    if (len == 0) return true;
    return std.mem.eql(u32, oracle.ptr.?[0..len], native);
}

fn enumSpanEqual(oracle: c_abi.U32Span, native: anytype) bool {
    const len: usize = @intCast(oracle.len);
    if (len != native.len) return false;
    if (len == 0) return true;
    const raw = oracle.ptr.?;
    for (native, 0..) |value, i| {
        if (raw[i] != @backingInt(value)) return false;
    }
    return true;
}

fn writeSummary(
    report: *std.Io.Writer,
    partition: corpus.Partition,
    tolerance: f64,
    counters: *const [public_observables.len]Counter,
    totals: *const Totals,
) !void {
    try report.print(
        "native-oracle-diff {t}: {d} members, T={d:.4} bond lengths\n",
        .{ partition, totals.members, tolerance },
    );
    if (totals.agreed_rejections != 0) {
        try report.print("  agreed rejections: {d}\n", .{totals.agreed_rejections});
    }
    if (totals.native_unsupported != 0) {
        try report.print("  native unsupported: {d}\n", .{totals.native_unsupported});
    }
    if (totals.status_mismatches != 0) {
        try report.print("  status mismatches: {d}\n", .{totals.status_mismatches});
    }
    for (public_observables, counters) |observable, counter| {
        try report.print(
            "  {t}: {d}/{d} matched",
            .{ observable, counter.compared - counter.mismatched, counter.compared },
        );
        if (counter.ceiling_applied != 0) {
            try report.print(", {d} at the enumerated ceiling", .{counter.ceiling_applied});
        }
        try report.print("\n", .{});
    }
    if (totals.max_deviation_member_len != 0) {
        try report.print(
            "  max coordinate deviation: {d:.6} bond lengths ({s})\n",
            .{ totals.max_deviation_bond_lengths, totals.worstMember() },
        );
    }
}

/// The machine-readable half. cgz-7v2.4.2 recalibrates T from this file and
/// cgz-7v2.4.7 takes its performance trigger from the same milestone, so the
/// numbers have to survive outside a build log.
fn writeBaseline(
    writer: *std.Io.Writer,
    partition: corpus.Partition,
    tolerance: f64,
    counters: *const [public_observables.len]Counter,
    totals: *const Totals,
) !void {
    try writer.print("# native-vs-oracle baseline, one line per observable\n", .{});
    try writer.print("# partition\t{t}\n", .{partition});
    try writer.print("# members\t{d}\n", .{totals.members});
    try writer.print("# tolerance_bond_lengths\t{d:.6}\n", .{tolerance});
    try writer.print(
        "# max_deviation_bond_lengths\t{d:.6}\t{s}\n",
        .{ totals.max_deviation_bond_lengths, totals.worstMember() },
    );
    try writer.print("observable\tcompared\tmismatched\tceiling_applied\n", .{});
    for (public_observables, counters) |observable, counter| {
        try writer.print("{t}\t{d}\t{d}\t{d}\n", .{
            observable,
            counter.compared,
            counter.mismatched,
            counter.ceiling_applied,
        });
    }
}

fn fatal(comptime format: []const u8, arguments: anytype) error{NativeOracleDiffFailed} {
    std.log.err(format, arguments);
    return error.NativeOracleDiffFailed;
}
