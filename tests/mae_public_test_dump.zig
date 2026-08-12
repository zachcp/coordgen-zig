//! Converts one pinned MAE fixture into the tiny, deterministic text protocol
//! consumed by conformance/public_fixture_test_rehost.cpp. This keeps the
//! public-test rehost independent of maeparser while using the same reviewed
//! Zig MAE reader as the fixture gate.

const std = @import("std");
const conformance = @import("conformance");

const max_fixture_bytes = 64 << 20;

pub fn main(init: std.process.Init) !void {
    const arena = init.arena.allocator();
    const gpa = init.gpa;
    const io = init.io;
    const args = try init.minimal.args.toSlice(arena);
    if (args.len != 3) return error.InvalidArguments;

    const source = try std.Io.Dir.cwd().readFileAlloc(
        io,
        args[1],
        gpa,
        .limited(max_fixture_bytes),
    );
    defer gpa.free(source);

    var diagnostic: conformance.mae.Diagnostic = .{};
    var reading = conformance.mae.parse(gpa, source, &diagnostic) catch |err| {
        std.log.err("{s}:{d}: {t} at '{s}'", .{
            args[1], diagnostic.line, err, diagnostic.token,
        });
        return err;
    };
    defer reading.deinit();

    var buffer: [4096]u8 = undefined;
    const output = try std.Io.Dir.cwd().createFile(io, args[2], .{});
    defer output.close(io);
    var file_writer = output.writer(io, &buffer);
    const writer = &file_writer.interface;

    for (reading.structures) |structure| {
        try writer.print("structure {d} {d}\n", .{
            structure.atoms.len,
            structure.bonds.len,
        });
        for (structure.atoms) |atom| {
            const x: f32 = @floatCast(atom.x);
            const y: f32 = @floatCast(atom.y);
            try writer.print("atom {d} {d} {d}\n", .{
                atom.atomic_number,
                @as(u32, @bitCast(x)),
                @as(u32, @bitCast(y)),
            });
        }
        for (structure.bonds) |bond| {
            try writer.print("bond {d} {d} {d}\n", .{
                bond.from,
                bond.to,
                bond.order,
            });
        }
        try writer.writeAll("end\n");
    }
    try writer.flush();
}
