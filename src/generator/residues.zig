const std = @import("std");
const core = @import("core");
const model = @import("model");

/// Allocator-owned translation of upstream's `groupResiduesInSSEs` result.
/// Groups and members retain chain-lexicographic/residue-number order.
pub const SecondaryStructures = struct {
    allocator: std.mem.Allocator,
    offsets: []u32,
    residues: []core.ids.ResidueId,

    pub fn deinit(self: *SecondaryStructures) void {
        self.allocator.free(self.residues);
        self.allocator.free(self.offsets);
        self.* = undefined;
    }

    pub fn count(self: SecondaryStructures) usize {
        return self.offsets.len - 1;
    }

    pub fn members(self: SecondaryStructures, index: usize) []const core.ids.ResidueId {
        std.debug.assert(index < self.count());
        return self.residues[self.offsets[index]..self.offsets[index + 1]];
    }
};

const ResidueOrderContext = struct {
    residues: []const model.Residue,
    chains: []const u8,
};

/// Divide residues into secondary-structure runs. Upstream first orders chains
/// lexicographically and each chain by residue number, then starts a new run
/// after a gap greater than three or for blank/empty chain identifiers.
pub fn groupSecondaryStructures(
    allocator: std.mem.Allocator,
    residues: []const model.Residue,
    chain_bytes: []const u8,
) core.errors.Error!SecondaryStructures {
    const ordered = allocator.alloc(core.ids.ResidueId, residues.len) catch return error.OutOfMemory;
    errdefer allocator.free(ordered);
    for (ordered, 0..) |*id, index| id.* = core.ids.ResidueId.fromIndex(@intCast(index));
    std.mem.sortUnstable(
        core.ids.ResidueId,
        ordered,
        ResidueOrderContext{ .residues = residues, .chains = chain_bytes },
        residueLessThan,
    );

    const all_offsets = allocator.alloc(u32, residues.len + 1) catch return error.OutOfMemory;
    defer allocator.free(all_offsets);
    var group_count: usize = 0;
    var previous: ?model.Residue = null;
    for (ordered, 0..) |id, ordered_index| {
        const residue = residues[id.index()];
        const chain = try residueChain(residue, chain_bytes);
        const starts_group = if (previous) |prior|
            !std.mem.eql(u8, try residueChain(prior, chain_bytes), chain) or
                residue.residue_number -| prior.residue_number > 3 or
                chain.len == 0 or std.mem.eql(u8, chain, " ")
        else
            true;
        if (starts_group) {
            all_offsets[group_count] = @intCast(ordered_index);
            group_count += 1;
        }
        previous = residue;
    }
    all_offsets[group_count] = @intCast(ordered.len);
    const offsets = allocator.dupe(u32, all_offsets[0 .. group_count + 1]) catch return error.OutOfMemory;
    return .{ .allocator = allocator, .offsets = offsets, .residues = ordered };
}

fn residueLessThan(
    context: ResidueOrderContext,
    first_id: core.ids.ResidueId,
    second_id: core.ids.ResidueId,
) bool {
    const first = context.residues[first_id.index()];
    const second = context.residues[second_id.index()];
    const first_chain = residueChain(first, context.chains) catch unreachable;
    const second_chain = residueChain(second, context.chains) catch unreachable;
    return switch (std.mem.order(u8, first_chain, second_chain)) {
        .lt => true,
        .gt => false,
        .eq => if (first.residue_number != second.residue_number)
            first.residue_number < second.residue_number
        else
            first.id.index() < second.id.index(),
    };
}

fn residueChain(residue: model.Residue, chain_bytes: []const u8) core.errors.Error![]const u8 {
    const start = residue.chain_start;
    const end = @as(usize, start) + residue.chain_len;
    if (end > chain_bytes.len) return error.InvalidMapping;
    return chain_bytes[start..end];
}

test "secondary structures follow chain order, residue gaps, and blank-chain splits" {
    const chains = "BBAAA  ";
    const input = [_]model.Residue{
        .{ .id = core.ids.ResidueId.fromIndex(0), .atom = core.ids.AtomId.fromIndex(0), .chain_start = 0, .chain_len = 1, .residue_number = 9 },
        .{ .id = core.ids.ResidueId.fromIndex(1), .atom = core.ids.AtomId.fromIndex(1), .chain_start = 1, .chain_len = 1, .residue_number = 2 },
        .{ .id = core.ids.ResidueId.fromIndex(2), .atom = core.ids.AtomId.fromIndex(2), .chain_start = 2, .chain_len = 1, .residue_number = 1 },
        .{ .id = core.ids.ResidueId.fromIndex(3), .atom = core.ids.AtomId.fromIndex(3), .chain_start = 3, .chain_len = 1, .residue_number = 3 },
        .{ .id = core.ids.ResidueId.fromIndex(4), .atom = core.ids.AtomId.fromIndex(4), .chain_start = 4, .chain_len = 1, .residue_number = 8 },
        .{ .id = core.ids.ResidueId.fromIndex(5), .atom = core.ids.AtomId.fromIndex(5), .chain_start = 5, .chain_len = 1, .residue_number = 1 },
        .{ .id = core.ids.ResidueId.fromIndex(6), .atom = core.ids.AtomId.fromIndex(6), .chain_start = 6, .chain_len = 1, .residue_number = 2 },
    };
    var groups = try groupSecondaryStructures(std.testing.allocator, &input, chains);
    defer groups.deinit();
    try std.testing.expectEqual(@as(usize, 6), groups.count());
    try std.testing.expectEqualSlices(core.ids.ResidueId, &.{core.ids.ResidueId.fromIndex(5)}, groups.members(0));
    try std.testing.expectEqualSlices(core.ids.ResidueId, &.{core.ids.ResidueId.fromIndex(6)}, groups.members(1));
    try std.testing.expectEqualSlices(core.ids.ResidueId, &.{
        core.ids.ResidueId.fromIndex(2),
        core.ids.ResidueId.fromIndex(3),
    }, groups.members(2));
    try std.testing.expectEqualSlices(core.ids.ResidueId, &.{core.ids.ResidueId.fromIndex(4)}, groups.members(3));
    try std.testing.expectEqualSlices(core.ids.ResidueId, &.{core.ids.ResidueId.fromIndex(1)}, groups.members(4));
    try std.testing.expectEqualSlices(core.ids.ResidueId, &.{core.ids.ResidueId.fromIndex(0)}, groups.members(5));
}
