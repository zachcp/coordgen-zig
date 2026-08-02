//! Prints a corpus partition in the canonical dump form, with no oracle
//! involved.
//!
//! The build diffs this against `conformance/smiles_reference.cpp`, which
//! prints the same partition as upstream's own SMILES parser builds it. That
//! is what keeps the committed drug-like tables honest.
//!
//! Usage: corpus-dump --partition NAME [--count N]

const std = @import("std");
const conformance = @import("conformance");

const corpus = conformance.corpus;

pub fn main(init: std.process.Init) !void {
    const arena = init.arena.allocator();
    const gpa = init.gpa;
    const io = init.io;

    var partition: ?corpus.Partition = null;
    var requested_count: u32 = 0;

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
        } else {
            return fatal("unknown option '{s}'", .{args[index]});
        }
    }

    const selected = partition orelse return fatal("--partition is required", .{});

    var buffer: [16 * 1024]u8 = undefined;
    var stdout_file: std.Io.File.Writer = .init(.stdout(), io, &buffer);
    const stdout = &stdout_file.interface;

    const count = selected.memberCount(requested_count);
    for (0..count) |raw_index| {
        const molecule = try corpus.generate(gpa, selected, @intCast(raw_index));
        defer molecule.deinit(gpa);
        try corpus.writeMolecule(stdout, molecule);
    }
    try stdout.flush();
}

fn fatal(comptime format: []const u8, arguments: anytype) error{CorpusDumpFailed} {
    std.log.err(format, arguments);
    return error.CorpusDumpFailed;
}
