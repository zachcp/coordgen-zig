const std = @import("std");

pub fn build(b: *std.Build) void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});
    const coordgen_dependency = b.dependency("coordgen", .{
        .target = target,
        .optimize = optimize,
    });

    // Use the example shipped by the dependency, not a copy in this consumer
    // package. This proves the example is present in the package artifact and
    // only reaches coordgen through its exported module.
    const consumer = b.createModule(.{
        .root_source_file = coordgen_dependency.path("examples/validate_input.zig"),
        .target = target,
        .optimize = optimize,
        .imports = &.{.{
            .name = "coordgen",
            .module = coordgen_dependency.module("coordgen"),
        }},
        .link_libc = false,
        .link_libcpp = false,
    });
    const example = b.addExecutable(.{
        .name = "coordgen-package-example",
        .root_module = consumer,
    });
    const run_example = b.addRunArtifact(example);
    run_example.expectExitCode(0);
    const test_step = b.step("test", "Run the example shipped by the coordgen package");
    test_step.dependOn(&run_example.step);
}
