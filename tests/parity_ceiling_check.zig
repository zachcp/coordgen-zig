//! Keeps `conformance/parity_ceiling.tsv` honest without an oracle.
//!
//! The ceiling enumeration is portable because it is keyed by corpus member
//! identity rather than by build, and `{partition}/{index}` is a portable
//! identity because `src/conformance/corpus.zig` generates every member from
//! nothing but those two values with an integer-only, value-pinned PRNG. That
//! argument holds across architecture, toolchain, and optimize mode; the one
//! thing it does not survive is an edit to the generator itself, which would
//! silently repoint every enumerated identity at a different molecule.
//!
//! So each row also carries the SHA-256 of the canonical dump of its member,
//! and this checker regenerates the member and compares. An identity whose
//! bytes moved fails by name instead of quietly excusing a molecule nobody
//! measured.
//!
//! Usage: parity-ceiling-check --ceiling FILE

const std = @import("std");
const conformance = @import("conformance");

const comparison = conformance.comparison;
const corpus = conformance.corpus;

const max_ceiling_bytes = 1 << 20;
/// The largest adversarial member is three components of at most thirty atoms
/// each, so its canonical dump is a few kilobytes.
const max_dump_bytes = 1 << 20;

pub fn main(init: std.process.Init) !void {
    const arena = init.arena.allocator();
    const gpa = init.gpa;
    const io = init.io;

    var ceiling_path: ?[]const u8 = null;
    const args = try init.minimal.args.toSlice(arena);
    var index: usize = 1;
    while (index < args.len) : (index += 2) {
        if (index + 1 >= args.len) return fatal("option '{s}' needs a value", .{args[index]});
        if (std.mem.eql(u8, args[index], "--ceiling")) {
            ceiling_path = args[index + 1];
        } else {
            return fatal("unknown option '{s}'", .{args[index]});
        }
    }
    const path = ceiling_path orelse return fatal("--ceiling is required", .{});

    const text = std.Io.Dir.cwd().readFileAlloc(io, path, arena, .limited(max_ceiling_bytes)) catch |err|
        return fatal("cannot read '{s}': {t}", .{ path, err });

    const table: comparison.CeilingTable = .{ .text = text };
    table.validateUnique() catch |err| return fatal("'{s}' is malformed: {t}", .{ path, err });

    const dump_buffer = try arena.alloc(u8, max_dump_bytes);

    var failures: usize = 0;
    var rows_seen: usize = 0;
    var digests: std.StringArrayHashMapUnmanaged([]const u8) = .empty;

    var rows = table.iterate();
    while (rows.next() catch |err| return fatal("'{s}' is malformed: {t}", .{ path, err })) |row| {
        rows_seen += 1;

        const split = std.mem.indexOfScalar(u8, row.member, '/').?;
        const partition_name = row.member[0..split];
        const partition = std.meta.stringToEnum(corpus.Partition, partition_name) orelse {
            std.log.err("member '{s}' names no corpus partition", .{row.member});
            failures += 1;
            continue;
        };
        const member_index = std.fmt.parseInt(u32, row.member[split + 1 ..], 10) catch {
            std.log.err("member '{s}' has no numeric index", .{row.member});
            failures += 1;
            continue;
        };

        // The digest belongs to the member, so every row naming that member
        // must agree. Disagreement means one of them was edited in isolation.
        const recorded = try digests.getOrPut(arena, row.member);
        if (recorded.found_existing) {
            if (!std.mem.eql(u8, recorded.value_ptr.*, row.input_sha256)) {
                std.log.err(
                    "member '{s}' carries two different input digests",
                    .{row.member},
                );
                failures += 1;
            }
            continue;
        }
        recorded.value_ptr.* = row.input_sha256;

        const molecule = corpus.generate(gpa, partition, member_index) catch |err| {
            std.log.err("member '{s}' does not generate: {t}", .{ row.member, err });
            failures += 1;
            continue;
        };
        defer molecule.deinit(gpa);

        var writer = std.Io.Writer.fixed(dump_buffer);
        try corpus.writeMolecule(&writer, molecule);

        var digest: [std.crypto.hash.sha2.Sha256.digest_length]u8 = undefined;
        std.crypto.hash.sha2.Sha256.hash(writer.buffered(), &digest, .{});
        var hexadecimal: [digest.len * 2]u8 = undefined;
        _ = std.fmt.bufPrint(&hexadecimal, "{x}", .{&digest}) catch unreachable;

        if (!std.mem.eql(u8, &hexadecimal, row.input_sha256)) {
            std.log.err(
                "member '{s}' generates {s}, but the ceiling records {s}; " ++
                    "the corpus generator changed and every enumerated identity " ++
                    "must be re-measured (cgz-r26)",
                .{ row.member, &hexadecimal, row.input_sha256 },
            );
            failures += 1;
        }
    }

    if (rows_seen == 0) return fatal("'{s}' enumerates nothing", .{path});
    if (failures != 0) return error.ParityCeilingCheckFailed;
}

fn fatal(comptime format: []const u8, arguments: anytype) error{ParityCeilingCheckFailed} {
    std.log.err(format, arguments);
    return error.ParityCeilingCheckFailed;
}
