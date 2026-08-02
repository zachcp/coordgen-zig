//! Runs a corpus partition through the pinned C++ oracle and dumps every
//! observable in a canonical, diffable form.
//!
//! One dump is produced per oracle build under comparison — native, a
//! cross-compiled architecture, and a build whose global `operator new` hands
//! out descending addresses. `tests/corpus_classify.zig` turns those dumps
//! into the stability classification.
//!
//! Order-sensitive collections are emitted twice: once in the order the
//! oracle produced them, and once canonically sorted. That distinction is the
//! point of the exercise — a permuted ring order and a different ring set are
//! different findings, and only the second one means the observable is
//! unusable for exact comparison.
//!
//! Usage: oracle-corpus-run --partition NAME --count N --output PATH

const std = @import("std");
const conformance = @import("conformance");
const c_abi = @import("c_abi");

const corpus = conformance.corpus;
const probe_types = conformance.probe_types;

extern fn coordgen_generate(input: *const c_abi.Input, result: *c_abi.Result) u32;
extern fn coordgen_result_free(result: *c_abi.Result) void;
extern fn coordgen_probe_generate(
    input: *const c_abi.Input,
    result: *probe_types.ProbeResult,
) u32;
extern fn coordgen_probe_result_free(result: *probe_types.ProbeResult) void;

/// Names match src/core/errors.zig's ErrorCode; the oracle adapter returns the
/// same values as the native library will.
fn errorName(code: u32) []const u8 {
    return switch (code) {
        0 => "ok",
        1 => "empty_graph",
        2 => "too_many_items",
        3 => "invalid_atomic_number",
        4 => "invalid_bond_order",
        5 => "invalid_atom_index",
        6 => "invalid_stereo",
        7 => "invalid_coordinate",
        8 => "invalid_option",
        9 => "invalid_mapping",
        10 => "out_of_memory",
        11 => "unsupported",
        12 => "internal",
        else => "unknown",
    };
}

pub fn main(init: std.process.Init) !void {
    const arena = init.arena.allocator();
    const gpa = init.gpa;
    const io = init.io;

    var partition: ?corpus.Partition = null;
    var requested_count: u32 = 0;
    var output_path: ?[]const u8 = null;

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
        } else if (std.mem.eql(u8, args[index], "--output")) {
            output_path = value;
        } else {
            return fatal("unknown option '{s}'", .{args[index]});
        }
    }

    const selected = partition orelse return fatal("--partition is required", .{});
    const path = output_path orelse return fatal("--output is required", .{});
    const count = selected.memberCount(requested_count);
    if (count == 0) return fatal("corpus partition '{t}' would be empty", .{selected});

    var file = try std.Io.Dir.cwd().createFile(io, path, .{});
    defer file.close(io);
    var buffer: [64 * 1024]u8 = undefined;
    var file_writer = file.writer(io, &buffer);
    const out = &file_writer.interface;

    try out.print("# corpus {t} members={d}\n", .{ selected, count });

    for (0..count) |raw_index| {
        const member_index: u32 = @intCast(raw_index);
        const molecule = try corpus.generate(gpa, selected, member_index);
        defer molecule.deinit(gpa);
        try emitMember(gpa, out, molecule);
        // Flushed per member so that a dump truncated by an oracle crash still
        // names the input that crashed it.
        try out.flush();
    }

    try out.flush();
}

fn emitMember(
    gpa: std.mem.Allocator,
    out: *std.Io.Writer,
    molecule: corpus.Molecule,
) !void {
    const key = MemberKey{ .partition = molecule.partition, .index = molecule.index };

    const atoms = try gpa.alloc(c_abi.AtomInput, molecule.atoms.len);
    defer gpa.free(atoms);
    for (atoms, molecule.atoms) |*atom, source| {
        atom.* = .{
            .atomic_number = source.atomic_number,
            .formal_charge = source.formal_charge,
        };
    }

    const bonds = try gpa.alloc(c_abi.BondInput, molecule.bonds.len);
    defer gpa.free(bonds);
    for (bonds, molecule.bonds) |*bond, source| {
        bond.* = .{ .start = source.start, .end = source.end, .order = source.order };
    }

    const input: c_abi.Input = .{
        .atoms = .{ .ptr = atoms.ptr, .len = @intCast(atoms.len) },
        .bonds = .{ .ptr = bonds.ptr, .len = @intCast(bonds.len) },
    };

    try out.print("input\t{f}\tatoms={d} bonds={d} bucket={t}\n", .{
        key,
        molecule.atoms.len,
        molecule.bonds.len,
        corpus.SizeBucket.of(molecule.atoms.len),
    });

    var result: c_abi.Result = .{};
    const status = coordgen_generate(&input, &result);
    defer if (status == 0) coordgen_result_free(&result);
    try out.print("status\t{f}\t{s}\n", .{ key, errorName(status) });
    if (status == 0) {
        try out.print("clean_pose\t{f}\t{d}\n", .{ key, result.clean_pose });
        try emitCoordinates(out, key, result);
        try emitU32Span(out, key, "input_to_internal", result.input_to_internal);
        try emitU32Span(out, key, "internal_to_input", result.internal_to_input);
        try emitU32Span(out, key, "effective_bond_orders", result.effective_bond_orders);
        try emitU32Span(out, key, "bond_displays", result.bond_displays);
        try emitU32Span(out, key, "atom_stereo", result.atom_stereo);
    }

    var probe: probe_types.ProbeResult = .{};
    const probe_status = coordgen_probe_generate(&input, &probe);
    defer if (probe_status == 0) coordgen_probe_result_free(&probe);
    try out.print("probe_status\t{f}\t{s}\n", .{ key, errorName(probe_status) });
    if (probe_status == 0) try emitProbe(gpa, out, key, probe);
}

const MemberKey = struct {
    partition: corpus.Partition,
    index: u32,

    pub fn format(self: MemberKey, writer: *std.Io.Writer) std.Io.Writer.Error!void {
        try writer.print("{t}/{d}", .{ self.partition, self.index });
    }
};

fn emitCoordinates(out: *std.Io.Writer, key: MemberKey, result: c_abi.Result) !void {
    // Raw bit patterns: a classification of "identical" must mean identical,
    // not "prints the same at six decimals".
    try out.print("coordinates\t{f}\t", .{key});
    const points = if (result.coordinates.ptr) |ptr| ptr[0..result.coordinates.len] else &.{};
    for (points, 0..) |point, position| {
        if (position != 0) try out.writeByte(' ');
        try out.print("{x:0>8},{x:0>8}", .{ @as(u32, @bitCast(point.x)), @as(u32, @bitCast(point.y)) });
    }
    try out.writeByte('\n');
}

fn emitU32Span(
    out: *std.Io.Writer,
    key: MemberKey,
    name: []const u8,
    span: c_abi.U32Span,
) !void {
    try out.print("{s}\t{f}\t", .{ name, key });
    const values = span.ptr orelse &[_]u32{};
    for (values[0..span.len], 0..) |value, position| {
        if (position != 0) try out.writeByte(' ');
        try out.print("{d}", .{value});
    }
    try out.writeByte('\n');
}

fn emitProbe(
    gpa: std.mem.Allocator,
    out: *std.Io.Writer,
    key: MemberKey,
    probe: probe_types.ProbeResult,
) !void {
    try out.print("probe_clean_pose\t{f}\t{d}\n", .{ key, probe.clean_pose });
    try emitU32Span(out, key, "morgan_ranks", probe.morgan_ranks);

    const ring_atoms = spanValues(probe.ring_atoms);
    const fragment_atoms = spanValues(probe.fragment_atoms);
    const fragment_rings = spanValues(probe.fragment_rings);
    const component_atoms = spanValues(probe.component_atoms);

    var rings: std.ArrayList([]const u8) = .empty;
    defer freeRecords(gpa, &rings);
    for (probe.ringSlice()) |ring| {
        try rings.append(gpa, try std.fmt.allocPrint(gpa, "[{f}]", .{
            U32List{ .values = slice(ring_atoms, ring.atom_start, ring.atom_count) },
        }));
    }
    try emitRecords(out, key, "rings", rings.items);

    var fragments: std.ArrayList([]const u8) = .empty;
    defer freeRecords(gpa, &fragments);
    for (probe.fragmentSlice()) |fragment| {
        try fragments.append(gpa, try std.fmt.allocPrint(
            gpa,
            "[parent={d} component={d} flags={d} template={d} atoms={f} rings={f}]",
            .{
                fragment.parent,
                fragment.component,
                fragment.flags,
                fragment.template_match,
                U32List{ .values = slice(fragment_atoms, fragment.atom_start, fragment.atom_count) },
                U32List{ .values = slice(fragment_rings, fragment.ring_start, fragment.ring_count) },
            },
        ));
    }
    try emitRecords(out, key, "fragments", fragments.items);

    // Structure and derived floats are separate observables. A penalty that
    // moves by one unit in the last place is a different finding from a DOF
    // that changed state, and lumping them together would report the first as
    // if it were the second.
    var dofs: std.ArrayList([]const u8) = .empty;
    defer freeRecords(gpa, &dofs);
    var dof_penalties: std.ArrayList([]const u8) = .empty;
    defer freeRecords(gpa, &dof_penalties);
    for (probe.dofSlice()) |dof| {
        try dofs.append(gpa, try std.fmt.allocPrint(
            gpa,
            "[kind={d} fragment={d} state={d}/{d}/{d} tier={d} atoms={d},{d} ring={d} multiplier={d}]",
            .{
                @backingInt(dof.kind),
                dof.fragment,
                dof.current_state,
                dof.optimal_state,
                dof.state_count,
                dof.tier,
                dof.atom_a,
                dof.atom_b,
                dof.ring,
                dof.variant_penalty_multiplier,
            },
        ));
        try dof_penalties.append(gpa, try std.fmt.allocPrint(gpa, "[{d}={x:0>8}]", .{
            dof.id,
            @as(u32, @bitCast(dof.current_penalty)),
        }));
    }
    try emitRecords(out, key, "dofs", dofs.items);
    try emitRecords(out, key, "dof_penalties", dof_penalties.items);

    var templates: std.ArrayList([]const u8) = .empty;
    defer freeRecords(gpa, &templates);
    for (probe.templateMappingSlice()) |mapping| {
        try templates.append(gpa, try std.fmt.allocPrint(gpa, "[{d}->{d}]", .{
            mapping.input_atom,
            mapping.template_atom,
        }));
    }
    try emitRecords(out, key, "template_mappings", templates.items);

    var components: std.ArrayList([]const u8) = .empty;
    defer freeRecords(gpa, &components);
    var component_transforms: std.ArrayList([]const u8) = .empty;
    defer freeRecords(gpa, &component_transforms);
    for (probe.componentSlice()) |component| {
        try components.append(gpa, try std.fmt.allocPrint(gpa, "[atoms={f} status={d}]", .{
            U32List{ .values = slice(component_atoms, component.atom_start, component.atom_count) },
            @backingInt(component.transform_status),
        }));
        try component_transforms.append(gpa, try std.fmt.allocPrint(
            gpa,
            "[{x:0>8},{x:0>8},{x:0>8},{x:0>8},{x:0>8},{x:0>8}]",
            .{
                @as(u32, @bitCast(component.transform[0])),
                @as(u32, @bitCast(component.transform[1])),
                @as(u32, @bitCast(component.transform[2])),
                @as(u32, @bitCast(component.transform[3])),
                @as(u32, @bitCast(component.transform[4])),
                @as(u32, @bitCast(component.transform[5])),
            },
        ));
    }
    try emitRecords(out, key, "components", components.items);
    try emitRecords(out, key, "component_transforms", component_transforms.items);
}

/// Emits a collection twice: `name` in oracle order, `name_set` sorted. A
/// difference in only the first means the contents survived but their order
/// did not.
fn emitRecords(
    out: *std.Io.Writer,
    key: MemberKey,
    name: []const u8,
    records: [][]const u8,
) !void {
    try out.print("{s}\t{f}\t", .{ name, key });
    for (records, 0..) |record, position| {
        if (position != 0) try out.writeByte(' ');
        try out.writeAll(record);
    }
    try out.writeByte('\n');

    std.mem.sort([]const u8, records, {}, lessThanBytes);
    try out.print("{s}_set\t{f}\t", .{ name, key });
    for (records, 0..) |record, position| {
        if (position != 0) try out.writeByte(' ');
        try out.writeAll(record);
    }
    try out.writeByte('\n');
}

fn lessThanBytes(_: void, left: []const u8, right: []const u8) bool {
    return std.mem.order(u8, left, right) == .lt;
}

fn freeRecords(gpa: std.mem.Allocator, records: *std.ArrayList([]const u8)) void {
    for (records.items) |record| gpa.free(record);
    records.deinit(gpa);
}

fn spanValues(span: c_abi.U32Span) []const u32 {
    const values = span.ptr orelse return &.{};
    return values[0..span.len];
}

fn slice(values: []const u32, start: u32, count: u32) []const u32 {
    if (start >= values.len) return &.{};
    const end = @min(values.len, start + count);
    return values[start..end];
}

const U32List = struct {
    values: []const u32,

    pub fn format(self: U32List, writer: *std.Io.Writer) std.Io.Writer.Error!void {
        for (self.values, 0..) |value, position| {
            if (position != 0) try writer.writeByte(',');
            try writer.print("{d}", .{value});
        }
    }
};

fn fatal(comptime format: []const u8, arguments: anytype) error{CorpusRunFailed} {
    std.log.err(format, arguments);
    return error.CorpusRunFailed;
}
