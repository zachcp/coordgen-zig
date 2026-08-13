const fragments = @import("layout/fragments.zig");
const basic = @import("layout/basic.zig");
pub const templates = @import("layout/templates.zig");
pub const macrocycle = @import("layout/macrocycle.zig");

pub const Fragmentation = fragments.Fragmentation;
pub const initializeCoordinates = basic.initializeCoordinates;
pub const initializeCoordinatesWithOptions = basic.initializeCoordinatesWithOptions;
pub const captureFragmentFrames = basic.captureFragmentFrames;

test {
    _ = fragments;
    _ = basic;
    _ = templates;
    _ = macrocycle;
}
