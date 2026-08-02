//! Loads the pinned upstream `.mae` fixtures with the minimal Zig reader and
//! checks each one against expectations recorded in the build graph.
//!
//! Usage: `mae-fixture-check (--fixture PATH --expect KEY=VALUE,...)...`
//!
//! Every expected field is a total that can be derived from the fixture text
//! by other tools, so the gate stays independent of this reader's own output.

const std = @import("std");
const conformance = @import("conformance");

const mae = conformance.mae;

const max_fixture_bytes = 64 << 20;

const Expectation = struct {
    structures: ?u32 = null,
    atoms: ?u32 = null,
    bonds: ?u32 = null,
    atomic_number_sum: ?u64 = null,
    bond_order_sum: ?u64 = null,
    x_sum: ?f64 = null,
    y_sum: ?f64 = null,
};

pub fn main(init: std.process.Init) !void {
    const arena = init.arena.allocator();
    const gpa = init.gpa;
    const io = init.io;

    var stdout_buffer: [4096]u8 = undefined;
    var stdout_file: std.Io.File.Writer = .init(.stdout(), io, &stdout_buffer);
    const stdout = &stdout_file.interface;

    const args = try init.minimal.args.toSlice(arena);
    if (args.len < 2) return fatal("usage: mae-fixture-check (--fixture PATH --expect KEY=VALUE,...)...", .{});

    var checked: usize = 0;
    var index: usize = 1;
    while (index < args.len) {
        if (!std.mem.eql(u8, args[index], "--fixture")) {
            return fatal("expected --fixture, found '{s}'", .{args[index]});
        }
        if (index + 3 > args.len) return fatal("--fixture needs a path and an --expect list", .{});
        const path = args[index + 1];
        if (!std.mem.eql(u8, args[index + 2], "--expect")) {
            return fatal("fixture '{s}' has no --expect list", .{path});
        }
        if (index + 4 > args.len) return fatal("--expect needs a KEY=VALUE list", .{});
        const expectation = try parseExpectation(path, args[index + 3]);
        index += 4;

        const source = std.Io.Dir.cwd().readFileAlloc(
            io,
            path,
            gpa,
            .limited(max_fixture_bytes),
        ) catch |err| return fatal("cannot read fixture '{s}': {t}", .{ path, err });
        defer gpa.free(source);

        var diagnostic: mae.Diagnostic = .{};
        var reading = mae.parse(gpa, source, &diagnostic) catch |err| return fatal(
            "{s}:{d}: {t} at '{s}'",
            .{ path, diagnostic.line, err, diagnostic.token },
        );
        defer reading.deinit();

        const summary = mae.summarize(reading.structures);
        try stdout.print(
            "{s}: structures={d} atoms={d} bonds={d} atomic_number_sum={d} " ++
                "bond_order_sum={d} x_sum={d:.6} y_sum={d:.6}\n",
            .{
                path,
                summary.structures,
                summary.atoms,
                summary.bonds,
                summary.atomic_number_sum,
                summary.bond_order_sum,
                summary.x_sum,
                summary.y_sum,
            },
        );
        try stdout.flush();
        try compare(path, expectation, summary);
        checked += 1;
    }

    try stdout.print("mae-fixture-check: {d} fixtures matched their recorded totals\n", .{checked});
    try stdout.flush();
}

fn parseExpectation(path: []const u8, list: []const u8) !Expectation {
    var expectation: Expectation = .{};
    var fields = std.mem.splitScalar(u8, list, ',');
    while (fields.next()) |field| {
        const split = std.mem.indexOfScalar(u8, field, '=') orelse
            return fatal("fixture '{s}': expectation '{s}' is not KEY=VALUE", .{ path, field });
        const key = field[0..split];
        const text = field[split + 1 ..];
        if (std.mem.eql(u8, key, "structures")) {
            expectation.structures = try unsigned(u32, path, field, text);
        } else if (std.mem.eql(u8, key, "atoms")) {
            expectation.atoms = try unsigned(u32, path, field, text);
        } else if (std.mem.eql(u8, key, "bonds")) {
            expectation.bonds = try unsigned(u32, path, field, text);
        } else if (std.mem.eql(u8, key, "atomic_number_sum")) {
            expectation.atomic_number_sum = try unsigned(u64, path, field, text);
        } else if (std.mem.eql(u8, key, "bond_order_sum")) {
            expectation.bond_order_sum = try unsigned(u64, path, field, text);
        } else if (std.mem.eql(u8, key, "x_sum")) {
            expectation.x_sum = try real(path, field, text);
        } else if (std.mem.eql(u8, key, "y_sum")) {
            expectation.y_sum = try real(path, field, text);
        } else {
            return fatal("fixture '{s}': unknown expectation key '{s}'", .{ path, key });
        }
    }
    return expectation;
}

fn unsigned(comptime T: type, path: []const u8, field: []const u8, text: []const u8) !T {
    return std.fmt.parseInt(T, text, 10) catch
        fatal("fixture '{s}': expectation '{s}' is not an integer", .{ path, field });
}

fn real(path: []const u8, field: []const u8, text: []const u8) !f64 {
    return std.fmt.parseFloat(f64, text) catch
        fatal("fixture '{s}': expectation '{s}' is not a real number", .{ path, field });
}

fn compare(path: []const u8, expectation: Expectation, summary: mae.Summary) !void {
    try expectCount(path, "structures", expectation.structures, summary.structures);
    try expectCount(path, "atoms", expectation.atoms, summary.atoms);
    try expectCount(path, "bonds", expectation.bonds, summary.bonds);
    try expectCount(path, "atomic_number_sum", expectation.atomic_number_sum, summary.atomic_number_sum);
    try expectCount(path, "bond_order_sum", expectation.bond_order_sum, summary.bond_order_sum);
    try expectReal(path, "x_sum", expectation.x_sum, summary.x_sum);
    try expectReal(path, "y_sum", expectation.y_sum, summary.y_sum);
}

fn expectCount(path: []const u8, name: []const u8, expected: anytype, actual: anytype) !void {
    const wanted: u64 = expected orelse return;
    if (wanted != actual) return fatal(
        "fixture '{s}': {s} is {d}, recorded expectation is {d}",
        .{ path, name, actual, wanted },
    );
}

fn expectReal(path: []const u8, name: []const u8, expected: ?f64, actual: f64) !void {
    const wanted = expected orelse return;
    // The fixtures carry six decimal places; this tolerance separates a
    // formatting difference from a value that actually changed.
    const tolerance = 5e-7 * @max(1.0, @abs(wanted));
    if (!(@abs(wanted - actual) <= tolerance)) return fatal(
        "fixture '{s}': {s} is {d:.6}, recorded expectation is {d:.6}",
        .{ path, name, actual, wanted },
    );
}

fn fatal(comptime format: []const u8, arguments: anytype) error{FixtureCheckFailed} {
    std.log.err(format, arguments);
    return error.FixtureCheckFailed;
}
