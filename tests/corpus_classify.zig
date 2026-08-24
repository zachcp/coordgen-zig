//! Turns oracle corpus dumps into the stability classification and checks it
//! against the recorded expectations.
//!
//! Three dumps of the same corpus are compared: an ascending-address baseline,
//! an ascending-address build for a different architecture, and a native build
//! whose global `operator new` hands out descending addresses. For every member
//! and every observable the classifier records whether that observable survived
//! each perturbation, which is the per-input-per-observable classification
//! cgz-r13 requires: coordinates going unstable must not silently demote ring
//! sets or fragment trees.
//!
//! What is gated and what is published are deliberately different things.
//! Coordinate values and their deviations are per-(architecture, toolchain,
//! optimize-mode) artifacts, so byte-comparing a committed manifest across
//! platforms would produce false failures. The gate is therefore the
//! portable claim - which observables must never move, and how often the rest
//! are allowed to - recorded in conformance/parity_expectations.tsv. The
//! manifest itself is published evidence, regenerated per build.
//!
//! The order axis carries a second, stricter gate. The fraction ceilings in
//! parity_expectations.tsv bound how often the oracle may disagree with itself
//! but never say on which members, so an aggregate under budget can hide a
//! newly unstable member. conformance/parity_ceiling.tsv enumerates the
//! permitted (member, observable) pairs instead, and this classifier fails on
//! any order-unstable pair missing from it (cgz-r13, cgz-r26).
//!
//! Usage: corpus-classify (--partition NAME --baseline FILE --architecture FILE
//!            --allocator-order FILE)... --expectations FILE --ceiling FILE
//!            [--manifest FILE]

const std = @import("std");
const conformance = @import("conformance");

const corpus = conformance.corpus;
const comparison = conformance.comparison;

/// One bond length in CoordGen units. Deviations are reported in these units
/// because a raw float distance means nothing without the scale.
const bond_length: f64 = 50;

const max_dump_bytes = 512 << 20;

const Axis = enum {
    architecture,
    allocator_order,

    fn label(self: Axis) []const u8 {
        return switch (self) {
            .architecture => "arch",
            .allocator_order => "order",
        };
    }
};

const Member = struct {
    atom_count: u32 = 0,
    bond_count: u32 = 0,
    bucket: []const u8 = "unknown",
    observables: std.StringArrayHashMapUnmanaged([]const u8) = .empty,
};

const Dump = struct {
    members: std.StringArrayHashMapUnmanaged(Member) = .empty,

    fn deinit(self: *Dump, allocator: std.mem.Allocator) void {
        for (self.members.values()) |*member| member.observables.deinit(allocator);
        self.members.deinit(allocator);
    }

    fn parse(arena: std.mem.Allocator, source: []const u8) !Dump {
        var dump: Dump = .{};
        var lines = std.mem.splitScalar(u8, source, '\n');
        while (lines.next()) |line| {
            if (line.len == 0 or line[0] == '#') continue;
            var fields = std.mem.splitScalar(u8, line, '\t');
            const observable = fields.next() orelse continue;
            const key = fields.next() orelse return error.MalformedDumpLine;
            const payload = fields.next() orelse return error.MalformedDumpLine;

            const entry = try dump.members.getOrPut(arena, key);
            if (!entry.found_existing) entry.value_ptr.* = .{};

            if (std.mem.eql(u8, observable, "input")) {
                entry.value_ptr.* = try parseInputRecord(payload, entry.value_ptr.*);
                continue;
            }
            try entry.value_ptr.observables.put(arena, observable, payload);
        }
        return dump;
    }

    fn parseInputRecord(payload: []const u8, current: Member) !Member {
        var member = current;
        var fields = std.mem.splitScalar(u8, payload, ' ');
        while (fields.next()) |field| {
            const split = std.mem.indexOfScalar(u8, field, '=') orelse continue;
            const key = field[0..split];
            const value = field[split + 1 ..];
            if (std.mem.eql(u8, key, "atoms")) {
                member.atom_count = try std.fmt.parseInt(u32, value, 10);
            } else if (std.mem.eql(u8, key, "bonds")) {
                member.bond_count = try std.fmt.parseInt(u32, value, 10);
            } else if (std.mem.eql(u8, key, "bucket")) {
                member.bucket = value;
            }
        }
        return member;
    }
};

/// Per (partition, observable) ceiling on how often an observable may fail to
/// survive an axis. An observable with no recorded expectation is an error,
/// not a pass: a new observable must be classified before it can be gated.
const Expectation = struct {
    partition: []const u8,
    observable: []const u8,
    max_unstable: [std.enums.values(Axis).len]f64,
};

const Classification = struct {
    key: []const u8,
    partition: []const u8,
    bucket: []const u8,
    atom_count: u32,
    status: []const u8,
    clean_pose: []const u8,
    unstable: [std.enums.values(Axis).len]std.ArrayList([]const u8),
    max_deviation: [std.enums.values(Axis).len]f64,
};

const Counter = struct {
    measured: u64 = 0,
    unstable: [std.enums.values(Axis).len]u64 = @splat(0),
};

/// The dump protocol is itself a conformance contract. Expectations place
/// ceilings on values; this list prevents an emitter regression from deleting
/// an observable before it reaches those ceilings.
const required_observables = [_][]const u8{
    "status",
    "clean_pose",
    "coordinates",
    "input_to_internal",
    "internal_to_input",
    "effective_bond_orders",
    "bond_displays",
    "atom_stereo",
    "probe_status",
    "probe_clean_pose",
    "morgan_ranks",
    "rings",
    "rings_set",
    "fragments",
    "fragments_set",
    "dofs",
    "dofs_set",
    "dof_penalties",
    "dof_penalties_set",
    "template_mappings",
    "template_mappings_set",
    "components",
    "components_set",
    "component_transforms",
    "component_transforms_set",
};

pub fn main(init: std.process.Init) !void {
    const arena = init.arena.allocator();
    const io = init.io;

    var inputs: std.ArrayList(PartitionInput) = .empty;
    var expectations_path: ?[]const u8 = null;
    var ceiling_path: ?[]const u8 = null;
    var manifest_path: ?[]const u8 = null;

    const args = try init.minimal.args.toSlice(arena);
    var index: usize = 1;
    while (index < args.len) : (index += 2) {
        if (index + 1 >= args.len) return fatal("option '{s}' needs a value", .{args[index]});
        const value = args[index + 1];
        if (std.mem.eql(u8, args[index], "--partition")) {
            try inputs.append(arena, .{ .name = value });
        } else if (std.mem.eql(u8, args[index], "--baseline")) {
            (try currentInput(inputs.items)).baseline = value;
        } else if (std.mem.eql(u8, args[index], "--architecture")) {
            (try currentInput(inputs.items)).architecture = value;
        } else if (std.mem.eql(u8, args[index], "--allocator-order")) {
            (try currentInput(inputs.items)).allocator_order = value;
        } else if (std.mem.eql(u8, args[index], "--expectations")) {
            expectations_path = value;
        } else if (std.mem.eql(u8, args[index], "--ceiling")) {
            ceiling_path = value;
        } else if (std.mem.eql(u8, args[index], "--manifest")) {
            manifest_path = value;
        } else {
            return fatal("unknown option '{s}'", .{args[index]});
        }
    }

    if (inputs.items.len == 0) return fatal("at least one --partition is required", .{});
    const expectations_file = expectations_path orelse
        return fatal("--expectations is required", .{});
    // Required, not optional: a gate that silently degrades to the fraction
    // ceilings when an argument is forgotten is the budget framing cgz-r13
    // rejected.
    const ceiling_file = ceiling_path orelse
        return fatal("--ceiling is required; the enumerated parity ceiling is not optional (cgz-r26)", .{});

    const expectations = try readExpectations(arena, io, expectations_file);
    const ceiling: comparison.CeilingTable = .{ .text = try readFile(arena, io, ceiling_file) };
    ceiling.validateUnique() catch |err|
        return fatal("ceiling file '{s}' is malformed: {t}", .{ ceiling_file, err });

    var classifications: std.ArrayList(Classification) = .empty;
    var counters: std.StringArrayHashMapUnmanaged(Counter) = .empty;

    for (inputs.items) |partition| {
        try classifyPartition(arena, io, partition, &classifications, &counters);
    }

    // Violations go to stderr, because a build step's captured stdout is not
    // shown when the step fails and a gate whose reason is invisible is
    // barely a gate. A passing run stays silent: the full per-observable
    // summary is written into the manifest either way.
    var report_buffer: [64 * 1024]u8 = undefined;
    var report_file: std.Io.File.Writer = .init(.stderr(), io, &report_buffer);
    const report = &report_file.interface;

    if (manifest_path) |path| {
        var file = try std.Io.Dir.cwd().createFile(io, path, .{});
        defer file.close(io);
        var buffer: [64 * 1024]u8 = undefined;
        var writer = file.writer(io, &buffer);
        try writeManifest(&writer.interface, classifications.items, counters);
        try writer.interface.flush();
    }

    const failures = try checkExpectations(report, counters, expectations) +
        try checkCeiling(report, classifications.items, ceiling);
    if (failures != 0) {
        try writeSummary(report, counters);
        try report.flush();
        return error.ParityExpectationsViolated;
    }
    try report.flush();
}

const PartitionInput = struct {
    name: []const u8,
    baseline: ?[]const u8 = null,
    architecture: ?[]const u8 = null,
    allocator_order: ?[]const u8 = null,
};

fn currentInput(items: []PartitionInput) !*PartitionInput {
    if (items.len == 0) return fatal("dump options must follow a --partition", .{});
    return &items[items.len - 1];
}

fn classifyPartition(
    arena: std.mem.Allocator,
    io: std.Io,
    partition: PartitionInput,
    classifications: *std.ArrayList(Classification),
    counters: *std.StringArrayHashMapUnmanaged(Counter),
) !void {
    const baseline_path = partition.baseline orelse
        return fatal("partition '{s}' has no --baseline dump", .{partition.name});
    const architecture_path = partition.architecture orelse
        return fatal("partition '{s}' has no --architecture dump", .{partition.name});
    const allocator_path = partition.allocator_order orelse
        return fatal("partition '{s}' has no --allocator-order dump", .{partition.name});

    const baseline = try Dump.parse(arena, try readFile(arena, io, baseline_path));
    var variants: [std.enums.values(Axis).len]Dump = undefined;
    variants[@backingInt(Axis.architecture)] =
        try Dump.parse(arena, try readFile(arena, io, architecture_path));
    variants[@backingInt(Axis.allocator_order)] =
        try Dump.parse(arena, try readFile(arena, io, allocator_path));

    if (baseline.members.count() == 0) {
        return fatal("partition '{s}' baseline dump has no members", .{partition.name});
    }
    for (std.enums.values(Axis)) |axis| {
        try validateDumpSchema(
            partition.name,
            baseline,
            variants[@backingInt(axis)],
            axis.label(),
            &required_observables,
            true,
        );
    }

    for (baseline.members.keys(), baseline.members.values()) |key, member| {
        var classification: Classification = .{
            .key = key,
            .partition = partition.name,
            .bucket = member.bucket,
            .atom_count = member.atom_count,
            .status = member.observables.get("status") orelse "missing",
            .clean_pose = member.observables.get("clean_pose") orelse "-",
            .unstable = .{ .empty, .empty },
            .max_deviation = @splat(0),
        };

        for (member.observables.keys(), member.observables.values()) |observable, payload| {
            const counter_key = try std.fmt.allocPrint(
                arena,
                "{s}\t{s}",
                .{ partition.name, observable },
            );
            const counter = try counters.getOrPut(arena, counter_key);
            if (!counter.found_existing) counter.value_ptr.* = .{};
            counter.value_ptr.measured += 1;

            for (std.enums.values(Axis)) |axis| {
                const slot = @backingInt(axis);
                const other = variants[slot].members.get(key) orelse {
                    return fatal(
                        "member '{s}' is missing from the {s} dump of partition '{s}'",
                        .{ key, axis.label(), partition.name },
                    );
                };
                const other_payload = other.observables.get(observable) orelse {
                    try classification.unstable[slot].append(arena, observable);
                    counter.value_ptr.unstable[slot] += 1;
                    continue;
                };
                if (std.mem.eql(u8, payload, other_payload)) continue;

                try classification.unstable[slot].append(arena, observable);
                counter.value_ptr.unstable[slot] += 1;
                if (std.mem.eql(u8, observable, "coordinates")) {
                    classification.max_deviation[slot] = try maxDeviation(payload, other_payload);
                }
            }
        }

        try classifications.append(arena, classification);
    }
}

fn validateDumpSchema(
    partition: []const u8,
    baseline: Dump,
    other: Dump,
    axis: []const u8,
    required: []const []const u8,
    log_failure: bool,
) !void {
    if (baseline.members.count() != other.members.count()) {
        return schemaFailure(
            log_failure,
            "partition '{s}' has {d} baseline members but {d} {s} members",
            .{ partition, baseline.members.count(), other.members.count(), axis },
        );
    }
    for (baseline.members.keys(), baseline.members.values()) |key, member| {
        const other_member = other.members.get(key) orelse return schemaFailure(
            log_failure,
            "member '{s}' is missing from the {s} dump of partition '{s}'",
            .{ key, axis, partition },
        );
        if (member.atom_count != other_member.atom_count or
            member.bond_count != other_member.bond_count or
            !std.mem.eql(u8, member.bucket, other_member.bucket))
        {
            return schemaFailure(
                log_failure,
                "member '{s}' input metadata differs in the {s} dump of partition '{s}'",
                .{ key, axis, partition },
            );
        }
        for (required) |observable| {
            if (!member.observables.contains(observable) or
                !other_member.observables.contains(observable))
            {
                return schemaFailure(
                    log_failure,
                    "member '{s}' is missing required observable '{s}' in the {s} dump of partition '{s}'",
                    .{ key, observable, axis, partition },
                );
            }
        }
        if (member.observables.count() != other_member.observables.count()) {
            return schemaFailure(
                log_failure,
                "member '{s}' has a different observable schema in the {s} dump of partition '{s}'",
                .{ key, axis, partition },
            );
        }
        for (member.observables.keys()) |observable| {
            if (!other_member.observables.contains(observable)) {
                return schemaFailure(
                    log_failure,
                    "member '{s}' is missing observable '{s}' in the {s} dump of partition '{s}'",
                    .{ key, observable, axis, partition },
                );
            }
        }
    }
}

fn schemaFailure(
    log_failure: bool,
    comptime format: []const u8,
    arguments: anytype,
) error{CorpusClassifyFailed} {
    if (log_failure) return fatal(format, arguments);
    return error.CorpusClassifyFailed;
}

/// Largest single-coordinate difference, in bond lengths. Both payloads are
/// raw f32 bit patterns, so this reports the size of a real divergence
/// without ever treating a formatting difference as one.
fn maxDeviation(left: []const u8, right: []const u8) !f64 {
    var left_points = std.mem.splitScalar(u8, left, ' ');
    var right_points = std.mem.splitScalar(u8, right, ' ');
    var largest: f64 = 0;
    while (left_points.next()) |left_point| {
        const right_point = right_points.next() orelse break;
        if (left_point.len == 0 or right_point.len == 0) continue;
        const left_split = std.mem.indexOfScalar(u8, left_point, ',') orelse return error.MalformedCoordinate;
        const right_split = std.mem.indexOfScalar(u8, right_point, ',') orelse return error.MalformedCoordinate;
        largest = @max(largest, try componentDeviation(left_point[0..left_split], right_point[0..right_split]));
        largest = @max(largest, try componentDeviation(left_point[left_split + 1 ..], right_point[right_split + 1 ..]));
    }
    return largest / bond_length;
}

fn componentDeviation(left: []const u8, right: []const u8) !f64 {
    const left_value: f32 = @bitCast(try std.fmt.parseInt(u32, left, 16));
    const right_value: f32 = @bitCast(try std.fmt.parseInt(u32, right, 16));
    if (std.math.isNan(left_value) and std.math.isNan(right_value)) return 0;
    return @abs(@as(f64, left_value) - @as(f64, right_value));
}

fn readFile(arena: std.mem.Allocator, io: std.Io, path: []const u8) ![]u8 {
    return std.Io.Dir.cwd().readFileAlloc(io, path, arena, .limited(max_dump_bytes)) catch |err|
        return fatal("cannot read '{s}': {t}", .{ path, err });
}

fn readExpectations(
    arena: std.mem.Allocator,
    io: std.Io,
    path: []const u8,
) ![]const Expectation {
    const source = try readFile(arena, io, path);
    var expectations: std.ArrayList(Expectation) = .empty;
    var lines = std.mem.splitScalar(u8, source, '\n');
    while (lines.next()) |line| {
        const trimmed = std.mem.trim(u8, line, " \t\r");
        if (trimmed.len == 0 or trimmed[0] == '#') continue;
        var fields = std.mem.splitScalar(u8, trimmed, '\t');
        const partition = fields.next() orelse return error.MalformedExpectation;
        const observable = fields.next() orelse return error.MalformedExpectation;
        const arch = fields.next() orelse return error.MalformedExpectation;
        const order = fields.next() orelse return error.MalformedExpectation;
        try expectations.append(arena, .{
            .partition = partition,
            .observable = observable,
            .max_unstable = .{
                std.fmt.parseFloat(f64, arch) catch return error.MalformedExpectation,
                std.fmt.parseFloat(f64, order) catch return error.MalformedExpectation,
            },
        });
    }
    if (expectations.items.len == 0) return fatal("expectations file '{s}' is empty", .{path});
    return expectations.items;
}

fn writeManifest(
    out: *std.Io.Writer,
    classifications: []const Classification,
    counters: std.StringArrayHashMapUnmanaged(Counter),
) !void {
    try out.writeAll(
        \\# Parity manifest: which comparison tier each corpus input and
        \\# observable is entitled to. Published evidence, regenerated by
        \\# `zig build parity-manifest -Denable-oracle=true`; the portable
        \\# gate lives in conformance/parity_expectations.tsv.
        \\#
        \\# Coordinate deviations are in bond lengths (one bond = 50 units)
        \\# and are per-(architecture, toolchain, optimize-mode) values: they
        \\# describe the build that produced this file, not a promise about
        \\# any other build.
        \\#
        \\
    );
    try out.writeAll("# member\tbucket\tatoms\tstatus\tclean_pose" ++
        "\tarch_unstable\torder_unstable\tarch_max_dev_bl\torder_max_dev_bl\n");
    for (classifications) |classification| {
        try out.print("{s}\t{s}\t{d}\t{s}\t{s}\t{f}\t{f}\t{d:.3}\t{d:.3}\n", .{
            classification.key,
            classification.bucket,
            classification.atom_count,
            classification.status,
            classification.clean_pose,
            ObservableList{ .items = classification.unstable[0].items },
            ObservableList{ .items = classification.unstable[1].items },
            classification.max_deviation[0],
            classification.max_deviation[1],
        });
    }

    try out.writeAll("#\n# partition\tobservable\tmeasured\tarch_unstable\torder_unstable\n");
    for (counters.keys(), counters.values()) |key, counter| {
        try out.print("# {s}\t{d}\t{d}\t{d}\n", .{
            key,
            counter.measured,
            counter.unstable[0],
            counter.unstable[1],
        });
    }
}

const ObservableList = struct {
    items: []const []const u8,

    pub fn format(self: ObservableList, writer: *std.Io.Writer) std.Io.Writer.Error!void {
        if (self.items.len == 0) return writer.writeAll("-");
        for (self.items, 0..) |item, position| {
            if (position != 0) try writer.writeByte(',');
            try writer.writeAll(item);
        }
    }
};

fn writeSummary(out: *std.Io.Writer, counters: std.StringArrayHashMapUnmanaged(Counter)) !void {
    try out.writeAll("corpus-classify: partition/observable measured arch-unstable order-unstable\n");
    for (counters.keys(), counters.values()) |key, counter| {
        try out.print("  {s}\t{d}\t{d}\t{d}\n", .{
            key,
            counter.measured,
            counter.unstable[0],
            counter.unstable[1],
        });
    }
}

fn checkExpectations(
    out: *std.Io.Writer,
    counters: std.StringArrayHashMapUnmanaged(Counter),
    expectations: []const Expectation,
) !usize {
    var failures: usize = 0;

    // An explicit expectation is a statement that an observable must be
    // measured. Without this check, a misspelled or accidentally removed
    // observable could leave a stale expectation row that no counter visits.
    for (expectations) |expectation| {
        if (std.mem.eql(u8, expectation.observable, "*")) continue;
        var key_buffer: [1024]u8 = undefined;
        const key = try std.fmt.bufPrint(
            &key_buffer,
            "{s}\t{s}",
            .{ expectation.partition, expectation.observable },
        );
        if (counters.contains(key)) continue;
        try out.print(
            "corpus-classify: recorded expectation {s}/{s} was not measured\n",
            .{ expectation.partition, expectation.observable },
        );
        failures += 1;
    }

    for (counters.keys(), counters.values()) |key, counter| {
        const split = std.mem.indexOfScalar(u8, key, '\t').?;
        const partition = key[0..split];
        const observable = key[split + 1 ..];

        const expectation = findExpectation(expectations, partition, observable) orelse {
            try out.print(
                "corpus-classify: {s}/{s} has no recorded expectation; classify it before gating on it\n",
                .{ partition, observable },
            );
            failures += 1;
            continue;
        };

        for (std.enums.values(Axis)) |axis| {
            const slot = @backingInt(axis);
            const measured: f64 = @floatFromInt(counter.measured);
            const fraction = @as(f64, @floatFromInt(counter.unstable[slot])) / measured;
            if (fraction <= expectation.max_unstable[slot]) continue;
            try out.print(
                "corpus-classify: {s}/{s} is {s}-unstable for {d}/{d} members ({d:.4}), above the recorded ceiling {d:.4}\n",
                .{
                    partition,
                    observable,
                    axis.label(),
                    counter.unstable[slot],
                    counter.measured,
                    fraction,
                    expectation.max_unstable[slot],
                },
            );
            failures += 1;
        }
    }
    return failures;
}

/// The enumerated half of the order-axis gate. `checkExpectations` bounds how
/// many pairs may be order-unstable; this bounds *which*, and the two run
/// together rather than one replacing the other.
///
/// Only the allocator-order axis is checked. The architecture axis is not a
/// parity ceiling: the differential runner executes oracle and native in one
/// process, build, and target, so architectural divergence is not a difference
/// it can observe (cgz-r13).
fn checkCeiling(
    out: *std.Io.Writer,
    classifications: []const Classification,
    ceiling: comparison.CeilingTable,
) !usize {
    var failures: usize = 0;
    const slot = @backingInt(Axis.allocator_order);

    for (classifications) |classification| {
        for (classification.unstable[slot].items) |observable_name| {
            const observable = std.meta.stringToEnum(comparison.Observable, observable_name) orelse {
                try out.print(
                    "corpus-classify: dumped observable '{s}' has no comparison.Observable; " ++
                        "the dump protocol and the comparison enum have drifted apart\n",
                    .{observable_name},
                );
                failures += 1;
                continue;
            };
            if (try ceiling.find(classification.key, observable) != null) continue;
            try out.print(
                "corpus-classify: {s}/{s} is order-unstable but is not enumerated in the " ++
                    "parity ceiling; the ceiling is an enumerated set of (member, observable) " ++
                    "pairs and never a budget, so this fails closed. Admitting it is a review " ++
                    "event: add the row with the measurement that placed it there " ++
                    "(cgz-r13, cgz-r26)\n",
                .{ classification.key, observable_name },
            );
            failures += 1;
        }
    }

    // A row is permission, not obligation, so an enumerated pair that happens
    // to be stable on this host is a pass. A row naming a member this corpus
    // does not contain is different: the file has stopped describing the
    // population it was measured on.
    var rows = ceiling.iterate();
    while (try rows.next()) |row| {
        var present = false;
        for (classifications) |classification| {
            if (std.mem.eql(u8, classification.key, row.member)) {
                present = true;
                break;
            }
        }
        if (present) continue;
        try out.print(
            "corpus-classify: the parity ceiling enumerates '{s}', which is not a member " ++
                "of the classified corpus\n",
            .{row.member},
        );
        failures += 1;
    }

    return failures;
}

fn findExpectation(
    expectations: []const Expectation,
    partition: []const u8,
    observable: []const u8,
) ?Expectation {
    var wildcard: ?Expectation = null;
    for (expectations) |expectation| {
        if (!std.mem.eql(u8, expectation.partition, partition)) continue;
        if (std.mem.eql(u8, expectation.observable, observable)) return expectation;
        if (std.mem.eql(u8, expectation.observable, "*")) wildcard = expectation;
    }
    return wildcard;
}

const test_digest = "0000000000000000000000000000000000000000000000000000000000000000";

test "an order-unstable pair outside the enumeration fails closed" {
    const classificationFor = struct {
        fn build(key: []const u8, order_unstable: std.ArrayList([]const u8)) Classification {
            return .{
                .key = key,
                .partition = "adversarial",
                .bucket = "medium",
                .atom_count = 36,
                .status = "ok",
                .clean_pose = "0",
                .unstable = .{ .empty, order_unstable },
                .max_deviation = @splat(0),
            };
        }
    }.build;

    const ceiling: comparison.CeilingTable = .{
        .text = "adversarial/1695\tcoordinates\t0.75\t" ++ test_digest ++ "\tcgz-r13\n",
    };

    var unstable_coordinates: std.ArrayList([]const u8) = .empty;
    defer unstable_coordinates.deinit(std.testing.allocator);
    try unstable_coordinates.append(std.testing.allocator, "coordinates");

    var buffer: [8192]u8 = undefined;
    var out = std.Io.Writer.fixed(&buffer);

    const enumerated = [_]Classification{classificationFor("adversarial/1695", unstable_coordinates)};
    try std.testing.expectEqual(@as(usize, 0), try checkCeiling(&out, &enumerated, ceiling));

    // The same observable on a member nobody measured is not excused by the
    // aggregate fraction having room in it.
    const other_member = [_]Classification{
        classificationFor("adversarial/1695", .empty),
        classificationFor("adversarial/7", unstable_coordinates),
    };
    try std.testing.expectEqual(@as(usize, 1), try checkCeiling(&out, &other_member, ceiling));

    // A row naming a member outside the corpus is a stale file, not a pass.
    const missing_member = [_]Classification{classificationFor("adversarial/7", .empty)};
    try std.testing.expectEqual(@as(usize, 1), try checkCeiling(&out, &missing_member, ceiling));
}

test "schema rejects observable, member, and input metadata drift" {
    const baseline_source =
        "input\tadversarial/0\tatoms=2 bonds=1 bucket=tiny\n" ++
        "status\tadversarial/0\tok\n";
    var baseline = try Dump.parse(std.testing.allocator, baseline_source);
    defer baseline.deinit(std.testing.allocator);

    var missing_observable = try Dump.parse(
        std.testing.allocator,
        "input\tadversarial/0\tatoms=2 bonds=1 bucket=tiny\n",
    );
    defer missing_observable.deinit(std.testing.allocator);
    try std.testing.expectError(
        error.CorpusClassifyFailed,
        validateDumpSchema("adversarial", baseline, missing_observable, "architecture", &.{"status"}, false),
    );

    var extra_member = try Dump.parse(
        std.testing.allocator,
        baseline_source ++
            "input\tadversarial/1\tatoms=2 bonds=1 bucket=tiny\n" ++
            "status\tadversarial/1\tok\n",
    );
    defer extra_member.deinit(std.testing.allocator);
    try std.testing.expectError(
        error.CorpusClassifyFailed,
        validateDumpSchema("adversarial", baseline, extra_member, "architecture", &.{"status"}, false),
    );

    const baseline_two_members = baseline_source ++
        "input\tadversarial/1\tatoms=2 bonds=1 bucket=tiny\n" ++
        "status\tadversarial/1\tok\n";
    var missing_member = try Dump.parse(std.testing.allocator, baseline_source);
    defer missing_member.deinit(std.testing.allocator);
    var two_members = try Dump.parse(std.testing.allocator, baseline_two_members);
    defer two_members.deinit(std.testing.allocator);
    try std.testing.expectError(
        error.CorpusClassifyFailed,
        validateDumpSchema("adversarial", two_members, missing_member, "architecture", &.{"status"}, false),
    );

    var changed_atoms = try Dump.parse(
        std.testing.allocator,
        "input\tadversarial/0\tatoms=3 bonds=1 bucket=tiny\n" ++
            "status\tadversarial/0\tok\n",
    );
    defer changed_atoms.deinit(std.testing.allocator);
    try std.testing.expectError(
        error.CorpusClassifyFailed,
        validateDumpSchema("adversarial", baseline, changed_atoms, "architecture", &.{"status"}, false),
    );

    var changed_bonds = try Dump.parse(
        std.testing.allocator,
        "input\tadversarial/0\tatoms=2 bonds=2 bucket=tiny\n" ++
            "status\tadversarial/0\tok\n",
    );
    defer changed_bonds.deinit(std.testing.allocator);
    try std.testing.expectError(
        error.CorpusClassifyFailed,
        validateDumpSchema("adversarial", baseline, changed_bonds, "architecture", &.{"status"}, false),
    );

    var changed_bucket = try Dump.parse(
        std.testing.allocator,
        "input\tadversarial/0\tatoms=2 bonds=1 bucket=small\n" ++
            "status\tadversarial/0\tok\n",
    );
    defer changed_bucket.deinit(std.testing.allocator);
    try std.testing.expectError(
        error.CorpusClassifyFailed,
        validateDumpSchema("adversarial", baseline, changed_bucket, "architecture", &.{"status"}, false),
    );
}

fn fatal(comptime format: []const u8, arguments: anytype) error{CorpusClassifyFailed} {
    std.log.err(format, arguments);
    return error.CorpusClassifyFailed;
}
