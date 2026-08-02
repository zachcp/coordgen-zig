const std = @import("std");

/// Public array positions are always fixed-width and remain in caller order.
pub const InputIndex = u32;
pub const invalid_input_index: InputIndex = std.math.maxInt(InputIndex);

fn Id(comptime name: []const u8) type {
    return enum(u32) {
        invalid = std.math.maxInt(u32),
        _,

        pub const kind_name = name;

        pub fn fromIndex(value: u32) @This() {
            std.debug.assert(value != std.math.maxInt(u32));
            return @fromBackingInt(@intCast(value));
        }

        pub fn index(self: @This()) u32 {
            return @backingInt(self);
        }

        pub fn isValid(self: @This()) bool {
            return self != .invalid;
        }
    };
}

/// Internal IDs are deliberately distinct enum(u32) types. They are never
/// interchangeable with public input indices without an explicit conversion.
pub const AtomId = Id("atom");
pub const BondId = Id("bond");
pub const RingId = Id("ring");
pub const FragmentId = Id("fragment");
pub const MoleculeId = Id("molecule");
pub const ResidueId = Id("residue");
pub const DofId = Id("dof");
pub const InteractionId = Id("interaction");

test "typed IDs have fixed representation and remain distinct" {
    try std.testing.expectEqual(@as(usize, 4), @sizeOf(AtomId));
    try std.testing.expectEqual(@as(usize, 4), @sizeOf(FragmentId));
    try std.testing.expectEqual(@as(u32, 7), AtomId.fromIndex(7).index());
    try std.testing.expect(AtomId != BondId);
}
