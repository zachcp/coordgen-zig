const fragments = @import("layout/fragments.zig");
const basic = @import("layout/basic.zig");

pub const Fragmentation = fragments.Fragmentation;
pub const initializeCoordinates = basic.initializeCoordinates;

test {
    _ = fragments;
    _ = basic;
}
