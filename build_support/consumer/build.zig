const std = @import("std");

pub fn build(b: *std.Build) void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});
    const coordgen_dependency = b.dependency("coordgen", .{
        .target = target,
        .optimize = optimize,
    });

    const consumer = b.createModule(.{
        .root_source_file = b.path("src/main.zig"),
        .target = target,
        .optimize = optimize,
        .imports = &.{.{
            .name = "coordgen",
            .module = coordgen_dependency.module("coordgen"),
        }},
        .link_libc = false,
        .link_libcpp = false,
    });
    const tests = b.addTest(.{ .root_module = consumer });
    const run_tests = b.addRunArtifact(tests);
    const test_step = b.step("test", "Compile and run an external package consumer");
    test_step.dependOn(&run_tests.step);
}
