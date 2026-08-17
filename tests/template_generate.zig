const std = @import("std");
const conformance = @import("conformance");
const fixture = @import("layout").templates.data;

const TemplateError = error{
    BondIndexOutOfRange,
    TemplateFieldOverflow,
};

fn sameCluster(dd: f32, representative: f32) bool {
    return @as(f64, dd) * 0.9 < @as(f64, representative) and
        @as(f64, dd) * 1.1 > @as(f64, representative);
}

fn normalize(atoms: []conformance.mae.Atom, bonds: []const conformance.mae.Bond) TemplateError!void {
    var dds: [256]f32 = undefined;
    var counts: [256]u32 = undefined;
    var cluster_count: usize = 0;
    for (bonds) |bond| {
        if (bond.from >= atoms.len or bond.to >= atoms.len) return error.BondIndexOutOfRange;
        const a = atoms[bond.from];
        const b = atoms[bond.to];
        const dx: f32 = @as(f32, @floatCast(a.x)) - @as(f32, @floatCast(b.x));
        const dy: f32 = @as(f32, @floatCast(a.y)) - @as(f32, @floatCast(b.y));
        const dd: f32 = dx * dx + dy * dy;
        var found = false;
        for (dds[0..cluster_count], 0..) |representative, index| {
            if (sameCluster(dd, representative)) {
                counts[index] += 1;
                found = true;
                break;
            }
        }
        if (!found) {
            if (cluster_count == dds.len) return error.TemplateFieldOverflow;
            dds[cluster_count] = dd;
            counts[cluster_count] = 1;
            cluster_count += 1;
        }
    }
    if (cluster_count == 0) return;
    var selected: usize = 0;
    for (counts[0..cluster_count], 0..) |count, index| if (count > counts[selected]) {
        selected = index;
    };
    const scale: f32 = @sqrt(dds[selected]);
    for (atoms) |*atom| {
        atom.x = @as(f32, @floatCast(atom.x)) / scale;
        atom.y = @as(f32, @floatCast(atom.y)) / scale;
    }
}

fn validateEncodable(structures: []const conformance.mae.Structure) TemplateError!void {
    var atom_total: usize = 0;
    var bond_total: usize = 0;
    for (structures) |structure| {
        if (structure.atoms.len > std.math.maxInt(u8) or structure.bonds.len > std.math.maxInt(u8)) {
            return error.TemplateFieldOverflow;
        }
        atom_total += structure.atoms.len;
        bond_total += structure.bonds.len;
        if (atom_total > std.math.maxInt(u16) or bond_total > std.math.maxInt(u16)) {
            return error.TemplateFieldOverflow;
        }
        for (structure.atoms) |atom| if (atom.atomic_number > std.math.maxInt(u8)) {
            return error.TemplateFieldOverflow;
        };
        for (structure.bonds) |bond| {
            if (bond.from >= structure.atoms.len or bond.to >= structure.atoms.len) return error.BondIndexOutOfRange;
            if (bond.from > std.math.maxInt(u8) or bond.to > std.math.maxInt(u8) or bond.order > std.math.maxInt(u8)) {
                return error.TemplateFieldOverflow;
            }
        }
    }
}

pub fn main(init: std.process.Init) !void {
    const io = init.io;
    const gpa = init.gpa;
    const args = try init.minimal.args.toSlice(gpa);
    defer gpa.free(args);
    if (args.len != 3 and args.len != 4) return error.InvalidArguments;
    const raw = args.len == 4 and std.mem.eql(u8, args[3], "--raw");
    if (args.len == 4 and !raw) return error.InvalidArguments;
    const source = try std.Io.Dir.cwd().readFileAlloc(io, args[1], gpa, .limited(16 * 1024 * 1024));
    defer gpa.free(source);
    var reading = try conformance.mae.parse(gpa, source, null);
    defer reading.deinit();
    if (reading.structures.len != 82) return error.UnexpectedTemplateCount;
    try validateEncodable(reading.structures);

    var file = try std.Io.Dir.cwd().createFile(io, args[2], .{});
    defer file.close(io);
    var buffer: [64 * 1024]u8 = undefined;
    var fw = file.writer(io, &buffer);
    const out = &fw.interface;
    try out.writeAll("// Generated template reference: do not edit.\n" ++
        "pub const Atom = struct { atomic_number: u8, x_bits: u32, y_bits: u32 };\n" ++
        "pub const Bond = struct { from: u8, to: u8, order: u8 };\n" ++
        "pub const Template = struct { atom_start: u16, atom_len: u8, bond_start: u16, bond_len: u8 };\n" ++
        "pub const atoms = [_]Atom{\n");
    var atom_total: usize = 0;
    var bond_total: usize = 0;
    for (reading.structures) |structure| {
        const atoms = @constCast(structure.atoms);
        if (!raw) try normalize(atoms, structure.bonds);
        for (atoms) |atom| try out.print("    .{{ .atomic_number = {d}, .x_bits = 0x{x:0>8}, .y_bits = 0x{x:0>8} }},\n", .{ atom.atomic_number, @as(u32, @bitCast(@as(f32, @floatCast(atom.x)))), @as(u32, @bitCast(@as(f32, @floatCast(atom.y)))) });
        atom_total += atoms.len;
        bond_total += structure.bonds.len;
    }
    if (atom_total != 1704 or bond_total != 1963) return error.UnexpectedTemplateTotals;
    try out.writeAll("};\npub const bonds = [_]Bond{\n");
    for (reading.structures) |structure| for (structure.bonds) |bond| try out.print("    .{{ .from = {d}, .to = {d}, .order = {d} }},\n", .{ bond.from, bond.to, bond.order });
    try out.writeAll("};\npub const templates = [_]Template{\n");
    var atom_start: usize = 0;
    var bond_start: usize = 0;
    for (reading.structures) |structure| {
        try out.print("    .{{ .atom_start = {d}, .atom_len = {d}, .bond_start = {d}, .bond_len = {d} }},\n", .{ atom_start, structure.atoms.len, bond_start, structure.bonds.len });
        atom_start += structure.atoms.len;
        bond_start += structure.bonds.len;
    }
    try out.writeAll("};\n");
    try out.flush();
}

test "normalization uses strict boundaries and first matching cluster" {
    var atoms = [_]conformance.mae.Atom{
        .{ .atomic_number = 6, .x = 0, .y = 0 },                    .{ .atomic_number = 6, .x = 1, .y = 0 },
        .{ .atomic_number = 6, .x = @sqrt(@as(f64, 0.9)), .y = 0 },
    };
    const bonds = [_]conformance.mae.Bond{ .{ .from = 0, .to = 1, .order = 1 }, .{ .from = 0, .to = 2, .order = 1 } };
    try normalize(&atoms, &bonds);
    try std.testing.expectEqual(@as(f32, 1), @as(f32, @floatCast(atoms[1].x)));
}

test "clustering promotes f32 values for upstream double-literal thresholds" {
    const dd: f32 = @bitCast(@as(u32, 0x3f0e38e9));
    const representative: f32 = @bitCast(@as(u32, 0x3f000005));
    try std.testing.expect(sameCluster(dd, representative));
    try std.testing.expect(!(dd * @as(f32, 0.9) < representative and dd * @as(f32, 1.1) > representative));
}

test "first cluster wins a count tie and empty bonds are unchanged" {
    var atoms = [_]conformance.mae.Atom{ .{ .atomic_number = 6, .x = 2, .y = 0 }, .{ .atomic_number = 6, .x = 4, .y = 0 }, .{ .atomic_number = 6, .x = 7, .y = 0 } };
    const bonds = [_]conformance.mae.Bond{ .{ .from = 0, .to = 1, .order = 1 }, .{ .from = 1, .to = 2, .order = 1 } };
    try normalize(&atoms, &bonds);
    try std.testing.expectEqual(@as(f32, 1), @as(f32, @floatCast(atoms[0].x)));
    var lone = [_]conformance.mae.Atom{.{ .atomic_number = 6, .x = 3, .y = 4 }};
    try normalize(&lone, &.{});
    try std.testing.expectEqual(@as(f64, 3), lone[0].x);
}

test "malformed bond index is rejected" {
    var atoms = [_]conformance.mae.Atom{.{ .atomic_number = 6, .x = 0, .y = 0 }};
    try std.testing.expectError(error.BondIndexOutOfRange, normalize(&atoms, &.{.{ .from = 0, .to = 1, .order = 1 }}));
}

test "committed fixture is valid immutable Zig data" {
    try std.testing.expectEqual(@as(usize, 82), fixture.templates.len);
    try std.testing.expectEqual(@as(usize, 1704), fixture.atoms.len);
    try std.testing.expectEqual(@as(usize, 1963), fixture.bonds.len);
    var atom_end: usize = 0;
    var bond_end: usize = 0;
    for (fixture.templates) |template| {
        try std.testing.expectEqual(atom_end, template.atom_start);
        try std.testing.expectEqual(bond_end, template.bond_start);
        atom_end += template.atom_len;
        bond_end += template.bond_len;
        for (fixture.bonds[template.bond_start..bond_end]) |bond| {
            try std.testing.expect(bond.from < template.atom_len);
            try std.testing.expect(bond.to < template.atom_len);
        }
    }
    try std.testing.expectEqual(fixture.atoms.len, atom_end);
    try std.testing.expectEqual(fixture.bonds.len, bond_end);
}

const TemplateSummary = struct {
    digest: u64,
    atom_count: usize,
    bond_count: usize,
};

fn mixTemplateField(digest: *u64, value: u64) void {
    digest.* = (digest.* ^ value) *% 0x100000001b3;
}

fn summarizeTemplates() TemplateSummary {
    var digest: u64 = 0xcbf29ce484222325;
    for (fixture.templates) |template| {
        mixTemplateField(&digest, template.atom_start);
        mixTemplateField(&digest, template.atom_len);
        mixTemplateField(&digest, template.bond_start);
        mixTemplateField(&digest, template.bond_len);
    }
    for (fixture.atoms) |atom| {
        mixTemplateField(&digest, atom.atomic_number);
        mixTemplateField(&digest, atom.x_bits);
        mixTemplateField(&digest, atom.y_bits);
    }
    for (fixture.bonds) |bond| {
        mixTemplateField(&digest, bond.from);
        mixTemplateField(&digest, bond.to);
        mixTemplateField(&digest, bond.order);
    }
    return .{
        .digest = digest,
        .atom_count = fixture.atoms.len,
        .bond_count = fixture.bonds.len,
    };
}

fn summarizeTemplatesOnThread(output: *TemplateSummary) void {
    output.* = summarizeTemplates();
}

test "immutable templates produce serial-identical results on two threads" {
    const serial = summarizeTemplates();
    var first: TemplateSummary = undefined;
    var second: TemplateSummary = undefined;
    const first_thread = try std.Thread.spawn(.{}, summarizeTemplatesOnThread, .{&first});
    const second_thread = try std.Thread.spawn(.{}, summarizeTemplatesOnThread, .{&second});
    first_thread.join();
    second_thread.join();

    try std.testing.expectEqual(serial, first);
    try std.testing.expectEqual(serial, second);
}
