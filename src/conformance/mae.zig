//! Minimal Maestro (`.mae`) reader for conformance fixtures.
//!
//! Upstream's test harness reaches the same data through maeparser and
//! `sketcherMaeReading.h`, which consumes exactly six fields:
//! `i_m_atomic_number`, `r_m_x_coord`, `r_m_y_coord` from `m_atom`, and
//! `i_m_from`, `i_m_to`, `i_m_order` from `m_bond`. Re-hosting that much on
//! the Zig side is what removes Boost and maeparser from the oracle build, so
//! this reader is deliberately scoped to those fields: every other block and
//! property is validated structurally and then discarded.
//!
//! This is not a general Maestro parser and is never part of an installed
//! artifact.

const std = @import("std");

/// Maestro block names this reader assigns meaning to. Everything else is
/// skipped after its structure is checked.
const structure_block = "f_m_ct";
const atom_block = "m_atom";
const bond_block = "m_bond";

pub const Atom = struct {
    atomic_number: u32,
    x: f64,
    y: f64,
};

pub const Bond = struct {
    /// Zero-based positions into `Structure.atoms`. Maestro stores atom
    /// references 1-based; the conversion happens once, here.
    from: u32,
    to: u32,
    order: u32,
};

pub const Structure = struct {
    atoms: []const Atom,
    bonds: []const Bond,
};

/// Owns every array it exposes through one arena. Exactly one `deinit`.
pub const Reading = struct {
    arena: std.heap.ArenaAllocator,
    structures: []const Structure,

    pub fn deinit(self: *Reading) void {
        self.arena.deinit();
        self.* = undefined;
    }
};

/// Optional parse-failure location. Populated only on error.
pub const Diagnostic = struct {
    line: usize = 0,
    token: []const u8 = "",
};

pub const ParseError = error{
    UnexpectedEndOfInput,
    UnexpectedToken,
    UnterminatedString,
    MissingBlockName,
    InvalidBlockSize,
    DuplicateProperty,
    MissingRequiredProperty,
    InvalidRowIndex,
    InvalidInteger,
    InvalidReal,
    MissingValue,
    DuplicateAtomBlock,
    DuplicateBondBlock,
    BlockOutsideStructure,
    NestedStructure,
    BondBeforeAtoms,
    AtomIndexOutOfRange,
} || std.mem.Allocator.Error;

pub fn parse(
    gpa: std.mem.Allocator,
    source: []const u8,
    diagnostic: ?*Diagnostic,
) ParseError!Reading {
    var arena: std.heap.ArenaAllocator = .init(gpa);
    errdefer arena.deinit();

    var parser: Parser = .{
        .arena = arena.allocator(),
        .tokenizer = .{ .source = source },
        .diagnostic = diagnostic,
    };
    try parser.parseDocument();

    return .{
        .arena = arena,
        .structures = try parser.structures.toOwnedSlice(parser.arena),
    };
}

/// Order-independent totals over a parsed file. Every field is derivable from
/// the fixture text by other tools, which is what makes it usable as a
/// conformance expectation rather than a self-confirming digest.
pub const Summary = struct {
    structures: u32 = 0,
    atoms: u32 = 0,
    bonds: u32 = 0,
    atomic_number_sum: u64 = 0,
    bond_order_sum: u64 = 0,
    x_sum: f64 = 0,
    y_sum: f64 = 0,
};

pub fn summarize(structures: []const Structure) Summary {
    var summary: Summary = .{ .structures = @intCast(structures.len) };
    for (structures) |structure| {
        summary.atoms += @intCast(structure.atoms.len);
        summary.bonds += @intCast(structure.bonds.len);
        for (structure.atoms) |atom| {
            summary.atomic_number_sum += atom.atomic_number;
            summary.x_sum += atom.x;
            summary.y_sum += atom.y;
        }
        for (structure.bonds) |bond| summary.bond_order_sum += bond.order;
    }
    return summary;
}

const Token = struct {
    kind: Kind,
    text: []const u8,
    /// Quoted values keep their escapes unresolved: no field this reader
    /// consumes is a string, so resolving them would be dead code that still
    /// had to be correct.
    quoted: bool = false,
    line: usize,

    const Kind = enum { value, open_brace, close_brace, separator, end };
};

const Tokenizer = struct {
    source: []const u8,
    index: usize = 0,
    line: usize = 1,

    fn next(self: *Tokenizer) ParseError!Token {
        self.skipTrivia();
        if (self.index == self.source.len) {
            return .{ .kind = .end, .text = "", .line = self.line };
        }

        const start_line = self.line;
        switch (self.source[self.index]) {
            '{' => {
                self.index += 1;
                return .{ .kind = .open_brace, .text = "{", .line = start_line };
            },
            '}' => {
                self.index += 1;
                return .{ .kind = .close_brace, .text = "}", .line = start_line };
            },
            '"' => return self.quotedValue(),
            else => {},
        }

        const start = self.index;
        while (self.index < self.source.len) : (self.index += 1) {
            switch (self.source[self.index]) {
                ' ', '\t', '\r', '\n', '{', '}', '"', '#' => break,
                else => {},
            }
        }
        const text = self.source[start..self.index];
        const kind: Token.Kind = if (std.mem.eql(u8, text, ":::")) .separator else .value;
        return .{ .kind = kind, .text = text, .line = start_line };
    }

    fn skipTrivia(self: *Tokenizer) void {
        while (self.index < self.source.len) {
            switch (self.source[self.index]) {
                '\n' => {
                    self.line += 1;
                    self.index += 1;
                },
                ' ', '\t', '\r' => self.index += 1,
                // Maestro comments run to end of line outside quoted values.
                '#' => while (self.index < self.source.len and
                    self.source[self.index] != '\n') : (self.index += 1)
                {},
                else => return,
            }
        }
    }

    fn quotedValue(self: *Tokenizer) ParseError!Token {
        const start_line = self.line;
        self.index += 1;
        const start = self.index;
        while (self.index < self.source.len) : (self.index += 1) {
            switch (self.source[self.index]) {
                '"' => {
                    const text = self.source[start..self.index];
                    self.index += 1;
                    return .{
                        .kind = .value,
                        .text = text,
                        .quoted = true,
                        .line = start_line,
                    };
                },
                '\\' => {
                    // Skip the escaped byte so an escaped quote does not end
                    // the value.
                    if (self.index + 1 == self.source.len) return error.UnterminatedString;
                    if (self.source[self.index + 1] == '\n') self.line += 1;
                    self.index += 1;
                },
                '\n' => self.line += 1,
                else => {},
            }
        }
        return error.UnterminatedString;
    }
};

const BlockName = struct {
    base: []const u8,
    /// Present for indexed blocks written as `name[count]`.
    count: ?u32,
};

const Builder = struct {
    atoms: std.ArrayList(Atom) = .empty,
    bonds: std.ArrayList(Bond) = .empty,
    seen_atom_block: bool = false,
    seen_bond_block: bool = false,
};

const Parser = struct {
    arena: std.mem.Allocator,
    tokenizer: Tokenizer,
    structures: std.ArrayList(Structure) = .empty,
    diagnostic: ?*Diagnostic,

    fn parseDocument(self: *Parser) ParseError!void {
        while (true) {
            // End of input is only legal here, between top-level blocks.
            const token = try self.tokenizer.next();
            switch (token.kind) {
                .end => return,
                // The `s_m_m2io_version` header block carries no name.
                .open_brace => try self.parseBlock("", null),
                .value => {
                    try self.expect(.open_brace);
                    if (std.mem.eql(u8, token.text, structure_block)) {
                        var builder: Builder = .{};
                        try self.parseBlock(token.text, &builder);
                        try self.structures.append(self.arena, .{
                            .atoms = try builder.atoms.toOwnedSlice(self.arena),
                            .bonds = try builder.bonds.toOwnedSlice(self.arena),
                        });
                    } else {
                        try self.parseBlock(token.text, null);
                    }
                },
                else => return self.fail(error.UnexpectedToken, token),
            }
        }
    }

    /// Parses one block body; the opening brace is already consumed.
    fn parseBlock(self: *Parser, raw_name: []const u8, builder: ?*Builder) ParseError!void {
        const name = try self.blockName(raw_name);
        const properties = try self.parsePropertyNames();

        if (name.count) |count| {
            try self.parseIndexedRows(name.base, count, properties, builder);
        } else {
            for (properties) |_| {
                const token = try self.next();
                if (token.kind != .value) return self.fail(error.UnexpectedToken, token);
            }
        }

        while (true) {
            const token = try self.next();
            switch (token.kind) {
                .close_brace => return,
                .value => {
                    try self.expect(.open_brace);
                    if (std.mem.eql(u8, token.text, structure_block)) {
                        return self.fail(error.NestedStructure, token);
                    }
                    try self.parseBlock(token.text, builder);
                },
                else => return self.fail(error.UnexpectedToken, token),
            }
        }
    }

    fn blockName(self: *Parser, raw_name: []const u8) ParseError!BlockName {
        const open = std.mem.indexOfScalar(u8, raw_name, '[') orelse
            return .{ .base = raw_name, .count = null };
        if (raw_name.len == 0 or raw_name[raw_name.len - 1] != ']') {
            return self.failText(error.InvalidBlockSize, raw_name);
        }
        const digits = raw_name[open + 1 .. raw_name.len - 1];
        const count = std.fmt.parseInt(u32, digits, 10) catch
            return self.failText(error.InvalidBlockSize, raw_name);
        if (open == 0) return self.failText(error.MissingBlockName, raw_name);
        return .{ .base = raw_name[0..open], .count = count };
    }

    fn parsePropertyNames(self: *Parser) ParseError![]const []const u8 {
        var names: std.ArrayList([]const u8) = .empty;
        while (true) {
            const token = try self.next();
            switch (token.kind) {
                .separator => break,
                .value => {
                    for (names.items) |existing| {
                        if (std.mem.eql(u8, existing, token.text)) {
                            return self.fail(error.DuplicateProperty, token);
                        }
                    }
                    try names.append(self.arena, token.text);
                },
                else => return self.fail(error.UnexpectedToken, token),
            }
        }
        return names.items;
    }

    fn parseIndexedRows(
        self: *Parser,
        base: []const u8,
        count: u32,
        properties: []const []const u8,
        builder: ?*Builder,
    ) ParseError!void {
        const is_atoms = std.mem.eql(u8, base, atom_block);
        const is_bonds = std.mem.eql(u8, base, bond_block);
        if ((is_atoms or is_bonds) and builder == null) {
            return self.failText(error.BlockOutsideStructure, base);
        }
        if (is_atoms and builder.?.seen_atom_block) {
            return self.failText(error.DuplicateAtomBlock, base);
        }
        if (is_bonds and builder.?.seen_bond_block) {
            return self.failText(error.DuplicateBondBlock, base);
        }
        if (is_bonds and !builder.?.seen_atom_block) {
            return self.failText(error.BondBeforeAtoms, base);
        }

        var columns: [3]usize = undefined;
        if (is_atoms) {
            columns = try self.requireColumns(properties, &.{
                "i_m_atomic_number",
                "r_m_x_coord",
                "r_m_y_coord",
            });
            builder.?.seen_atom_block = true;
            try builder.?.atoms.ensureUnusedCapacity(self.arena, count);
        } else if (is_bonds) {
            columns = try self.requireColumns(properties, &.{
                "i_m_from",
                "i_m_to",
                "i_m_order",
            });
            builder.?.seen_bond_block = true;
            try builder.?.bonds.ensureUnusedCapacity(self.arena, count);
        }

        // One row is the leading Maestro index plus one value per property.
        const row = try self.arena.alloc(Token, properties.len + 1);
        defer self.arena.free(row);

        for (0..count) |ordinal| {
            for (row) |*cell| {
                const token = try self.next();
                if (token.kind != .value) return self.fail(error.UnexpectedToken, token);
                cell.* = token;
            }
            const index = try self.integer(row[0]);
            if (index != ordinal + 1) return self.fail(error.InvalidRowIndex, row[0]);

            if (is_atoms) {
                builder.?.atoms.appendAssumeCapacity(.{
                    .atomic_number = try self.integer(row[columns[0] + 1]),
                    .x = try self.real(row[columns[1] + 1]),
                    .y = try self.real(row[columns[2] + 1]),
                });
            } else if (is_bonds) {
                const atom_count: u32 = @intCast(builder.?.atoms.items.len);
                const from = try self.atomReference(row[columns[0] + 1], atom_count);
                const to = try self.atomReference(row[columns[1] + 1], atom_count);
                builder.?.bonds.appendAssumeCapacity(.{
                    .from = from,
                    .to = to,
                    .order = try self.integer(row[columns[2] + 1]),
                });
            }
        }

        try self.expect(.separator);
    }

    fn requireColumns(
        self: *Parser,
        properties: []const []const u8,
        required: []const []const u8,
    ) ParseError![3]usize {
        std.debug.assert(required.len == 3);
        var columns: [3]usize = undefined;
        for (required, 0..) |name, slot| {
            columns[slot] = for (properties, 0..) |property, column| {
                if (std.mem.eql(u8, property, name)) break column;
            } else return self.failText(error.MissingRequiredProperty, name);
        }
        return columns;
    }

    fn atomReference(self: *Parser, token: Token, atom_count: u32) ParseError!u32 {
        const value = try self.integer(token);
        if (value == 0 or value > atom_count) {
            return self.fail(error.AtomIndexOutOfRange, token);
        }
        return value - 1;
    }

    fn integer(self: *Parser, token: Token) ParseError!u32 {
        if (isMissing(token)) return self.fail(error.MissingValue, token);
        return std.fmt.parseInt(u32, token.text, 10) catch
            self.fail(error.InvalidInteger, token);
    }

    fn real(self: *Parser, token: Token) ParseError!f64 {
        if (isMissing(token)) return self.fail(error.MissingValue, token);
        return std.fmt.parseFloat(f64, token.text) catch
            self.fail(error.InvalidReal, token);
    }

    fn isMissing(token: Token) bool {
        return !token.quoted and std.mem.eql(u8, token.text, "<>");
    }

    fn next(self: *Parser) ParseError!Token {
        const token = try self.tokenizer.next();
        if (token.kind == .end) return self.fail(error.UnexpectedEndOfInput, token);
        return token;
    }

    fn expect(self: *Parser, kind: Token.Kind) ParseError!void {
        const token = try self.next();
        if (token.kind != kind) return self.fail(error.UnexpectedToken, token);
    }

    fn fail(self: *Parser, err: ParseError, token: Token) ParseError {
        if (self.diagnostic) |diagnostic| {
            diagnostic.* = .{ .line = token.line, .token = token.text };
        }
        return err;
    }

    fn failText(self: *Parser, err: ParseError, text: []const u8) ParseError {
        if (self.diagnostic) |diagnostic| {
            diagnostic.* = .{ .line = self.tokenizer.line, .token = text };
        }
        return err;
    }
};

const testing = std.testing;

const two_atom_source =
    \\{
    \\ s_m_m2io_version
    \\ :::
    \\ 2.0.0
    \\}
    \\
    \\f_m_ct {
    \\ s_m_title
    \\ :::
    \\ ""
    \\ m_atom[2] {
    \\  # First column is atom index #
    \\  i_m_mmod_type
    \\  r_m_x_coord
    \\  r_m_y_coord
    \\  r_m_z_coord
    \\  i_m_atomic_number
    \\  :::
    \\  1 25 -1.500000 0.250000 0.000000 7
    \\  2 3 2.000000 -0.750000 0.000000 6
    \\  :::
    \\ }
    \\ m_bond[1] {
    \\  # First column is bond index #
    \\  i_m_from
    \\  i_m_to
    \\  i_m_order
    \\  :::
    \\  1 1 2 1
    \\  :::
    \\ }
    \\}
    \\
;

test "reads the six fields upstream consumes and rebases atom references" {
    var reading = try parse(testing.allocator, two_atom_source, null);
    defer reading.deinit();

    try testing.expectEqual(@as(usize, 1), reading.structures.len);
    const structure = reading.structures[0];
    try testing.expectEqual(@as(usize, 2), structure.atoms.len);
    try testing.expectEqual(@as(u32, 7), structure.atoms[0].atomic_number);
    try testing.expectEqual(@as(f64, -1.5), structure.atoms[0].x);
    try testing.expectEqual(@as(f64, 0.25), structure.atoms[0].y);
    try testing.expectEqual(@as(u32, 6), structure.atoms[1].atomic_number);

    try testing.expectEqual(@as(usize, 1), structure.bonds.len);
    try testing.expectEqual(Bond{ .from = 0, .to = 1, .order = 1 }, structure.bonds[0]);
}

test "summary totals are independent of structure order" {
    var reading = try parse(testing.allocator, two_atom_source, null);
    defer reading.deinit();

    const summary = summarize(reading.structures);
    try testing.expectEqual(@as(u32, 1), summary.structures);
    try testing.expectEqual(@as(u32, 2), summary.atoms);
    try testing.expectEqual(@as(u32, 1), summary.bonds);
    try testing.expectEqual(@as(u64, 13), summary.atomic_number_sum);
    try testing.expectEqual(@as(u64, 1), summary.bond_order_sum);
    try testing.expectEqual(@as(f64, 0.5), summary.x_sum);
    try testing.expectEqual(@as(f64, -0.5), summary.y_sum);
}

test "multiple structures and a structure without bonds" {
    const source =
        \\f_m_ct {
        \\ s_m_title
        \\ :::
        \\ "no bonds"
        \\ m_atom[1] {
        \\  i_m_atomic_number
        \\  r_m_x_coord
        \\  r_m_y_coord
        \\  :::
        \\  1 11 0.000000 0.000000
        \\  :::
        \\ }
        \\}
        \\f_m_ct {
        \\ s_m_title
        \\ :::
        \\ "second"
        \\ m_atom[2] {
        \\  i_m_atomic_number
        \\  r_m_x_coord
        \\  r_m_y_coord
        \\  :::
        \\  1 6 0.000000 0.000000
        \\  2 8 1.000000 0.000000
        \\  :::
        \\ }
        \\ m_bond[1] {
        \\  i_m_from
        \\  i_m_to
        \\  i_m_order
        \\  :::
        \\  1 2 1 2
        \\  :::
        \\ }
        \\}
    ;
    var reading = try parse(testing.allocator, source, null);
    defer reading.deinit();

    try testing.expectEqual(@as(usize, 2), reading.structures.len);
    try testing.expectEqual(@as(usize, 0), reading.structures[0].bonds.len);
    try testing.expectEqual(
        Bond{ .from = 1, .to = 0, .order = 2 },
        reading.structures[1].bonds[0],
    );
}

test "quoted values, escapes, and comments never end a block early" {
    const source =
        \\f_m_ct {
        \\ s_m_title
        \\ s_m_entry_name
        \\ :::
        \\ "a } brace and a \" quote"
        \\ "# not a comment"
        \\ # a real comment with a } brace
        \\ m_atom[1] {
        \\  i_m_atomic_number
        \\  r_m_x_coord
        \\  r_m_y_coord
        \\  s_m_label
        \\  :::
        \\  1 6 0.000000 0.000000 <>
        \\  :::
        \\ }
        \\}
    ;
    var reading = try parse(testing.allocator, source, null);
    defer reading.deinit();

    try testing.expectEqual(@as(usize, 1), reading.structures.len);
    try testing.expectEqual(@as(u32, 6), reading.structures[0].atoms[0].atomic_number);
}

test "malformed input is rejected with a line diagnostic" {
    const cases = [_]struct { source: []const u8, expected: anyerror }{
        .{ .source =
        \\f_m_ct {
        \\ :::
        \\ m_atom[1] {
        \\  i_m_atomic_number
        \\  r_m_x_coord
        \\  :::
        \\  1 6 0.000000
        \\  :::
        \\ }
        \\}
        , .expected = error.MissingRequiredProperty },
        .{ .source =
        \\f_m_ct {
        \\ :::
        \\ m_atom[1] {
        \\  i_m_atomic_number
        \\  r_m_x_coord
        \\  r_m_y_coord
        \\  :::
        \\  1 6 0.000000 0.000000
        \\  :::
        \\ }
        \\ m_bond[1] {
        \\  i_m_from
        \\  i_m_to
        \\  i_m_order
        \\  :::
        \\  1 1 3 1
        \\  :::
        \\ }
        \\}
        , .expected = error.AtomIndexOutOfRange },
        .{ .source =
        \\f_m_ct {
        \\ :::
        \\ m_atom[2] {
        \\  i_m_atomic_number
        \\  r_m_x_coord
        \\  r_m_y_coord
        \\  :::
        \\  1 6 0.000000 0.000000
        \\  3 6 0.000000 0.000000
        \\  :::
        \\ }
        \\}
        , .expected = error.InvalidRowIndex },
        .{ .source =
        \\f_m_ct {
        \\ :::
        \\ m_atom[1] {
        \\  i_m_atomic_number
        \\  r_m_x_coord
        \\  r_m_y_coord
        \\  :::
        \\  1 <> 0.000000 0.000000
        \\  :::
        \\ }
        \\}
        , .expected = error.MissingValue },
        .{ .source =
        \\f_m_ct {
        \\ :::
        \\ m_atom[1] {
        \\  i_m_atomic_number
        \\  r_m_x_coord
        \\  r_m_y_coord
        \\  :::
        \\  1 6 0.000000 0.000000
        \\  :::
        , .expected = error.UnexpectedEndOfInput },
        .{ .source =
        \\f_m_ct {
        \\ :::
        \\ f_m_ct {
        \\  :::
        \\ }
        \\}
        , .expected = error.NestedStructure },
        .{ .source =
        \\f_m_ct {
        \\ :::
        \\ m_bond[1] {
        \\  i_m_from
        \\  i_m_to
        \\  i_m_order
        \\  :::
        \\  1 1 2 1
        \\  :::
        \\ }
        \\}
        , .expected = error.BondBeforeAtoms },
    };

    for (cases) |case| {
        var diagnostic: Diagnostic = .{};
        try testing.expectError(case.expected, parse(testing.allocator, case.source, &diagnostic));
        try testing.expect(diagnostic.line > 0);
    }
}

test "a quoted value spanning physical lines does not desynchronize the reader" {
    // test/macrocycle.mae stores a multi-line quoted property value whose
    // continuation lines begin with digits at column 0, so a line-oriented
    // reader would mistake them for data rows.
    const source =
        \\f_m_ct {
        \\ s_sd_PUBCHEM\_CACTVS\_SUBSKEYS_(undefined)
        \\ r_sd_Quick_Properties_(MW)
        \\ :::
        \\ "18  19  3
        \\20  21  4
        \\48  56  8"
        \\ 342.4
        \\ m_atom[1] {
        \\  i_m_atomic_number
        \\  r_m_x_coord
        \\  r_m_y_coord
        \\  :::
        \\  1 6 -0.000000 0.000000
        \\  :::
        \\ }
        \\}
    ;
    var reading = try parse(testing.allocator, source, null);
    defer reading.deinit();

    try testing.expectEqual(@as(usize, 1), reading.structures.len);
    try testing.expectEqual(@as(usize, 1), reading.structures[0].atoms.len);
    try testing.expectEqual(@as(f64, 0), reading.structures[0].atoms[0].x);
}

fn parseAndDiscard(allocator: std.mem.Allocator, source: []const u8) !void {
    var reading = try parse(allocator, source, null);
    reading.deinit();
}

test "parser reports out of memory instead of truncating a fixture" {
    try testing.checkAllAllocationFailures(
        testing.allocator,
        parseAndDiscard,
        .{two_atom_source},
    );
}
