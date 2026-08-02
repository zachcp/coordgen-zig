const std = @import("std");
const layer_contract = @import("src/module_layers.zig");

fn createInternalModule(
    b: *std.Build,
    target: std.Build.ResolvedTarget,
    optimize: std.builtin.OptimizeMode,
    root_source_file: ?[]const u8,
) *std.Build.Module {
    return b.createModule(.{
        .root_source_file = if (root_source_file) |path| b.path(path) else null,
        .target = target,
        .optimize = optimize,
        .link_libc = false,
        .link_libcpp = false,
    });
}

fn wireApprovedModuleEdges(modules: []const *std.Build.Module) void {
    layer_contract.validate() catch |err| {
        std.debug.panic("invalid approved module graph: {s}", .{@errorName(err)});
    };

    const layer_count = @typeInfo(layer_contract.Layer).@"enum".field_names.len;
    if (modules.len != layer_count) @panic("build module table does not cover every declared layer");

    // The build edge set is not a second hand-maintained list: every addImport
    // comes directly from the reviewed allow-list in src/module_layers.zig.
    for (layer_contract.approved_edges) |edge| {
        const from = modules[@backingInt(edge.from)];
        const to = modules[@backingInt(edge.to)];
        from.addImport(@tagName(edge.to), to);
    }
}

pub fn build(b: *std.Build) !void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});

    const core = createInternalModule(b, target, optimize, "src/core.zig");
    const model = createInternalModule(b, target, optimize, "src/model.zig");
    const geometry = createInternalModule(b, target, optimize, "src/geometry.zig");
    const topology = createInternalModule(b, target, optimize, "build_support/empty_module.zig");
    const layout = createInternalModule(b, target, optimize, "build_support/empty_module.zig");
    const optimize_module = createInternalModule(b, target, optimize, "build_support/empty_module.zig");
    const generator = createInternalModule(b, target, optimize, "build_support/empty_module.zig");
    const api = createInternalModule(b, target, optimize, "src/api.zig");
    const c_abi = createInternalModule(b, target, optimize, "src/c_abi_types.zig");
    const conformance = createInternalModule(b, target, optimize, "src/conformance/probe_types.zig");
    const layered_modules = [_]*std.Build.Module{
        core,
        model,
        geometry,
        topology,
        layout,
        optimize_module,
        generator,
        api,
        c_abi,
        conformance,
    };
    wireApprovedModuleEdges(&layered_modules);

    const module_layers = createInternalModule(b, target, optimize, "src/module_layers.zig");

    // Production is deliberately std-only. Explicit false values make an
    // accidental libc/libc++ dependency a build error rather than a silent
    // expansion of the runtime contract.
    const coordgen = b.addModule("coordgen", .{
        .root_source_file = b.path("src/coordgen.zig"),
        .target = target,
        .optimize = optimize,
        .imports = &.{
            .{ .name = "core", .module = core },
            .{ .name = "model", .module = model },
            .{ .name = "geometry", .module = geometry },
            .{ .name = "api", .module = api },
            .{ .name = "c_abi", .module = c_abi },
            .{ .name = "module_layers", .module = module_layers },
        },
        .link_libc = false,
        .link_libcpp = false,
    });

    const library = b.addLibrary(.{
        .name = "coordgen",
        .root_module = coordgen,
    });
    library.installHeader(b.path("include/coordgen_abi.h"), "coordgen_abi.h");
    b.installArtifact(library);

    const module_tests = b.addTest(.{ .root_module = coordgen });
    const run_module_tests = b.addRunArtifact(module_tests);
    const conformance_tests = b.addTest(.{
        .name = "conformance-types-test",
        .root_module = conformance,
    });
    const run_conformance_tests = b.addRunArtifact(conformance_tests);
    const layer_tests = b.addTest(.{
        .name = "module-layer-test",
        .root_module = module_layers,
    });
    const run_layer_tests = b.addRunArtifact(layer_tests);

    const consumer_module = b.createModule(.{
        .root_source_file = b.path("build_support/package_consumer.zig"),
        .target = target,
        .optimize = optimize,
        .imports = &.{.{ .name = "coordgen", .module = coordgen }},
        .link_libc = false,
        .link_libcpp = false,
    });
    const consumer_tests = b.addTest(.{
        .name = "package-consumer-test",
        .root_module = consumer_module,
    });
    const run_consumer_tests = b.addRunArtifact(consumer_tests);

    const test_step = b.step("test", "Run native and package-consumer tests");
    test_step.dependOn(&run_module_tests.step);
    test_step.dependOn(&run_conformance_tests.step);
    test_step.dependOn(&run_layer_tests.step);
    test_step.dependOn(&run_consumer_tests.step);

    const module_graph_step = b.step("module-graph-check", "Validate the approved named-module edge set");
    module_graph_step.dependOn(&run_layer_tests.step);

    const policy_command = b.addSystemCommand(&.{ "sh", "tools/check-build-policy" });
    const external_consumer = b.addRunFile(std.Build.LazyPath.zig_exe);
    external_consumer.addArgs(&.{ "build", "test", "--summary", "failures", "--cache-dir" });
    external_consumer.addDirectoryArg(std.Build.LazyPath.cache_root.path(b, "external-consumer"));
    external_consumer.setCwd(b.path("build_support/consumer"));
    const package_step = b.step("package-check", "Verify package and dependency isolation policy");
    package_step.dependOn(&policy_command.step);
    package_step.dependOn(&run_consumer_tests.step);
    package_step.dependOn(&external_consumer.step);

    const abi_layout_module = b.createModule(.{
        .target = target,
        .optimize = optimize,
        .link_libc = true,
    });
    abi_layout_module.addIncludePath(b.path("include"));
    abi_layout_module.addCSourceFile(.{
        .file = b.path("tests/abi_layout.c"),
        .flags = &.{ "-std=c11", "-Wall", "-Wextra", "-Werror", "-pedantic" },
    });
    const abi_layout = b.addExecutable(.{
        .name = "abi-layout-test",
        .root_module = abi_layout_module,
    });
    const run_abi_layout = b.addRunArtifact(abi_layout);

    const abi_cpp_module = b.createModule(.{
        .target = target,
        .optimize = optimize,
        .link_libcpp = true,
    });
    abi_cpp_module.addIncludePath(b.path("include"));
    abi_cpp_module.addCSourceFile(.{
        .file = b.path("tests/abi_cpp_consumer.cpp"),
        .flags = &.{ "-std=c++17", "-Wall", "-Wextra", "-Werror" },
    });
    const abi_cpp_consumer = b.addObject(.{
        .name = "abi-cpp-consumer",
        .root_module = abi_cpp_module,
    });

    const probe_layout_module = b.createModule(.{
        .target = target,
        .optimize = optimize,
        .link_libc = true,
    });
    probe_layout_module.addIncludePath(b.path("conformance/include"));
    probe_layout_module.addCSourceFile(.{
        .file = b.path("tests/probe_layout.c"),
        .flags = &.{ "-std=c11", "-Wall", "-Wextra", "-Werror", "-pedantic" },
    });
    const probe_layout = b.addObject(.{
        .name = "probe-layout-test",
        .root_module = probe_layout_module,
    });

    const abi_step = b.step("abi-check", "Compile and run stable ABI and non-installed probe layout checks");
    abi_step.dependOn(&run_abi_layout.step);
    abi_step.dependOn(&abi_cpp_consumer.step);
    abi_step.dependOn(&probe_layout.step);
    abi_step.dependOn(&policy_command.step);

    // Maker configuration cannot inspect requested top-level step arguments.
    // The explicit option is therefore the stable opt-in boundary: without it,
    // dependencyLazy is never called and the oracle cannot be fetched.
    const enable_oracle = b.option(
        bool,
        "enable-oracle",
        "Fetch the pinned C++ oracle for conformance-only build steps",
    ) orelse false;

    const oracle_step = b.step("upstream-oracle", "Validate the pinned conformance-only oracle package");
    const conformance_step = b.step("conformance", "Run oracle/native conformance checks");

    if (enable_oracle) {
        const oracle = try b.dependencyLazy("coordgenlibs_oracle", .{});
        // URL archive fetching strips the single GitHub top-level directory;
        // the archive root itself is verified by tools/verify-upstream.
        const license_check = b.addCheckFile(oracle.path("LICENSE"), .{
            .expected_matches = &.{
                "BSD 3-Clause License",
                "Copyright (c) 2017, Schrödinger, Inc.",
            },
            .max_bytes = 4096,
        });
        oracle_step.dependOn(&license_check.step);
        conformance_step.dependOn(oracle_step);
    } else {
        const oracle_disabled = b.addFail(
            "oracle steps require -Denable-oracle=true; ordinary builds never fetch conformance sources",
        );
        oracle_step.dependOn(&oracle_disabled.step);
        conformance_step.dependOn(&oracle_disabled.step);
    }

    // Reserve the public step vocabulary without creating false-green gates.
    // Each owning bead removes its addFail dependency when it attaches the
    // corresponding artifacts.
    const examples_step = b.step("examples", "Build examples when example sources are present");
    const examples_pending = b.addFail("examples is not implemented yet; see cgz-7v2");
    examples_step.dependOn(&examples_pending.step);

    const regeneration_step = b.step("template-regeneration-check", "Check generated template data when available");
    const regeneration_pending = b.addFail("template-regeneration-check is not implemented yet; see cgz-r09");
    regeneration_step.dependOn(&regeneration_pending.step);

    const fuzz_step = b.step("fuzz", "Run the platform-selected fuzz harness when available");
    const fuzz_pending = b.addFail("fuzz is not implemented yet; see cgz-r15");
    fuzz_step.dependOn(&fuzz_pending.step);
}
