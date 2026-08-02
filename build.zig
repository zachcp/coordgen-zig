const std = @import("std");
const layer_contract = @import("src/module_layers.zig");

/// Every `.mae` file in the pinned oracle package: upstream's six test
/// fixtures plus the template source. Totals are derived from the fixture
/// text by an independent tokenizer, never from this repository's reader, so
/// the fixture gate cannot confirm itself.
const fixture_expectations = [_]struct { path: []const u8, totals: []const u8 }{
    .{
        .path = "test/test.mae",
        .totals = "structures=1,atoms=26,bonds=26,atomic_number_sum=128," ++
            "bond_order_sum=31,x_sum=0.000000,y_sum=0.000000",
    },
    .{
        .path = "test/test_mol.mae",
        .totals = "structures=1,atoms=27,bonds=31,atomic_number_sum=174," ++
            "bond_order_sum=40,x_sum=-3.763010,y_sum=-3.038242",
    },
    .{
        .path = "test/testChirality.mae",
        .totals = "structures=1,atoms=9,bonds=8,atomic_number_sum=50," ++
            "bond_order_sum=8,x_sum=3.348949,y_sum=0.213781",
    },
    .{
        .path = "test/macrocycle.mae",
        .totals = "structures=1,atoms=59,bonds=67,atomic_number_sum=384," ++
            "bond_order_sum=81,x_sum=0.000300,y_sum=-0.000300",
    },
    .{
        .path = "test/metalZobs.mae",
        .totals = "structures=1,atoms=6,bonds=5,atomic_number_sum=24," ++
            "bond_order_sum=5,x_sum=4.980544,y_sum=-0.475011",
    },
    .{
        .path = "test/nonterminalMetalZobs.mae",
        .totals = "structures=1,atoms=6,bonds=5,atomic_number_sum=24," ++
            "bond_order_sum=5,x_sum=4.980544,y_sum=-0.475011",
    },
    .{
        .path = "templates.mae",
        .totals = "structures=82,atoms=1704,bonds=1963,atomic_number_sum=10610," ++
            "bond_order_sum=2162,x_sum=-932.150248,y_sum=-316.790238",
    },
};

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
    const conformance = createInternalModule(b, target, optimize, "src/conformance.zig");
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
        .name = "conformance-module-test",
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
    probe_layout_module.addIncludePath(b.path("include"));
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

        // The upstream library is Boost-free and maeparser-free: only its
        // CMake test target needs them, which is why the harness and fixture
        // reading are re-hosted on the Zig side instead. Listing the pinned
        // translation units explicitly makes an upstream file appearing or
        // disappearing a build failure rather than a silent change in what
        // the oracle contains.
        const instrument_fragment_builder = b.addSystemCommand(&.{
            "sh",
            "conformance/patches/apply-instrumentation.sh",
        });
        instrument_fragment_builder.addFileArg(oracle.path("CoordgenFragmentBuilder.cpp"));
        instrument_fragment_builder.addFileArg(b.path("conformance/patches/CoordgenFragmentBuilder.instrumentation.patch"));
        const instrumented_fragment_builder = instrument_fragment_builder.addOutputFileArg("CoordgenFragmentBuilder.cpp");
        const instrument_minimizer = b.addSystemCommand(&.{
            "sh",
            "conformance/patches/apply-instrumentation.sh",
        });
        instrument_minimizer.addFileArg(oracle.path("sketcherMinimizer.cpp"));
        instrument_minimizer.addFileArg(b.path("conformance/patches/sketcherMinimizer.instrumentation.patch"));
        const instrumented_minimizer = instrument_minimizer.addOutputFileArg("sketcherMinimizer.cpp");

        const oracle_module = b.createModule(.{
            .target = target,
            .optimize = optimize,
            .link_libcpp = true,
        });
        oracle_module.addIncludePath(b.path("conformance/include"));
        oracle_module.addCSourceFiles(.{
            .root = oracle.path("."),
            .files = &.{
                "CoordgenFragmenter.cpp",
                "CoordgenMacrocycleBuilder.cpp",
                "CoordgenMinimizer.cpp",
                "CoordgenTemplates.cpp",
                "sketcherMinimizerAtom.cpp",
                "sketcherMinimizerBond.cpp",
                "sketcherMinimizerFragment.cpp",
                "sketcherMinimizerMarchingSquares.cpp",
                "sketcherMinimizerMolecule.cpp",
                "sketcherMinimizerResidue.cpp",
                "sketcherMinimizerResidueInteraction.cpp",
                "sketcherMinimizerRing.cpp",
            },
            .flags = &.{ "-std=c++17", "-DIN_COORDGEN" },
        });
        oracle_module.addCSourceFile(.{
            .file = instrumented_fragment_builder,
            .flags = &.{ "-std=c++17", "-DIN_COORDGEN" },
        });
        oracle_module.addCSourceFile(.{
            .file = instrumented_minimizer,
            .flags = &.{ "-std=c++17", "-DIN_COORDGEN" },
        });
        oracle_module.addCSourceFile(.{
            .file = b.path("conformance/oracle_hook.cpp"),
            .flags = &.{ "-std=c++17", "-DIN_COORDGEN" },
        });
        const oracle_library = b.addLibrary(.{
            .name = "coordgen-oracle",
            .linkage = .static,
            .root_module = oracle_module,
        });
        // Deliberately never installed: the oracle exists for conformance
        // comparison only.
        oracle_step.dependOn(&oracle_library.step);

        // The adapter has the production-shaped C ABI plus the deliberately
        // broader probe API.  It is conformance-only: neither artifact nor
        // header is installed with the native Zig library.
        const oracle_abi_module = b.createModule(.{
            .target = target,
            .optimize = optimize,
            .link_libcpp = true,
        });
        oracle_abi_module.addIncludePath(b.path("include"));
        oracle_abi_module.addIncludePath(b.path("conformance/include"));
        oracle_abi_module.addIncludePath(oracle.path("."));
        oracle_abi_module.addCSourceFile(.{
            .file = b.path("conformance/oracle_adapter.cpp"),
            .flags = &.{ "-std=c++17", "-Wall", "-Wextra", "-Werror" },
        });
        oracle_abi_module.linkLibrary(oracle_library);
        const oracle_abi = b.addLibrary(.{
            .name = "coordgen-oracle-abi",
            .linkage = .static,
            .root_module = oracle_abi_module,
        });
        oracle_step.dependOn(&oracle_abi.step);

        const oracle_abi_smoke_module = b.createModule(.{
            .target = target,
            .optimize = optimize,
            .link_libc = true,
            .link_libcpp = true,
        });
        oracle_abi_smoke_module.addIncludePath(b.path("include"));
        oracle_abi_smoke_module.addIncludePath(b.path("conformance/include"));
        oracle_abi_smoke_module.addCSourceFile(.{
            .file = b.path("tests/oracle_abi_smoke.c"),
            .flags = &.{ "-std=c11", "-Wall", "-Wextra", "-Werror", "-pedantic" },
        });
        oracle_abi_smoke_module.linkLibrary(oracle_abi);
        const oracle_abi_smoke = b.addExecutable(.{
            .name = "oracle-abi-smoke",
            .root_module = oracle_abi_smoke_module,
        });
        const run_oracle_abi_smoke = b.addRunArtifact(oracle_abi_smoke);
        run_oracle_abi_smoke.expectExitCode(0);
        oracle_step.dependOn(&run_oracle_abi_smoke.step);

        const oracle_abi_cpp_smoke_module = b.createModule(.{
            .target = target,
            .optimize = optimize,
            .link_libcpp = true,
        });
        oracle_abi_cpp_smoke_module.addIncludePath(b.path("include"));
        oracle_abi_cpp_smoke_module.addCSourceFile(.{
            .file = b.path("tests/oracle_abi_cpp_smoke.cpp"),
            .flags = &.{ "-std=c++17", "-Wall", "-Wextra", "-Werror" },
        });
        oracle_abi_cpp_smoke_module.linkLibrary(oracle_abi);
        const oracle_abi_cpp_smoke = b.addExecutable(.{
            .name = "oracle-abi-cpp-smoke",
            .root_module = oracle_abi_cpp_smoke_module,
        });
        const run_oracle_abi_cpp_smoke = b.addRunArtifact(oracle_abi_cpp_smoke);
        run_oracle_abi_cpp_smoke.expectExitCode(0);
        oracle_step.dependOn(&run_oracle_abi_cpp_smoke.step);

        // Upstream's own example program is the cheapest proof that the
        // static library links and runs. Its two atoms must come out one
        // bond-length (50) apart on the x axis, which is the baseline
        // recorded when this build was first verified against the pin.
        const oracle_example_module = b.createModule(.{
            .target = target,
            .optimize = optimize,
            .link_libcpp = true,
        });
        oracle_example_module.addCSourceFile(.{
            .file = oracle.path("example_dir/example.cpp"),
            .flags = &.{"-std=c++17"},
        });
        oracle_example_module.linkLibrary(oracle_library);
        const oracle_example = b.addExecutable(.{
            .name = "oracle-example",
            .root_module = oracle_example_module,
        });
        const run_oracle_example = b.addRunArtifact(oracle_example);
        run_oracle_example.expectExitCode(0);
        run_oracle_example.expectStdErrEqual("(-50, 0)  (0, 0)\n");
        oracle_step.dependOn(&run_oracle_example.step);

        // Fixture loading runs against the pinned fixtures themselves, so the
        // reader is exercised by the same bytes the oracle harness will use.
        const fixture_check_module = b.createModule(.{
            .root_source_file = b.path("tests/mae_fixture_check.zig"),
            .target = target,
            .optimize = optimize,
            .imports = &.{.{ .name = "conformance", .module = conformance }},
            .link_libc = false,
            .link_libcpp = false,
        });
        const fixture_check = b.addExecutable(.{
            .name = "mae-fixture-check",
            .root_module = fixture_check_module,
        });
        const run_fixture_check = b.addRunArtifact(fixture_check);
        run_fixture_check.expectExitCode(0);
        // Totals are derived from the fixture text independently of this
        // reader; see docs/architecture/CONFORMANCE_ORACLE.md.
        for (fixture_expectations) |fixture| {
            run_fixture_check.addArg("--fixture");
            run_fixture_check.addFileArg(oracle.path(fixture.path));
            run_fixture_check.addArg("--expect");
            run_fixture_check.addArg(fixture.totals);
        }
        oracle_step.dependOn(&run_fixture_check.step);

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
