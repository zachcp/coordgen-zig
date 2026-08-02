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

/// Corpus partitions and their pinned sizes. The recorded classification is
/// only meaningful for the exact population it was measured on, so the count
/// is part of the contract, not a tuning knob: 2000 adversarial members is
/// the size behind the divergence measurements in cgz-r05.
const corpus_partitions = [_]struct { name: []const u8, count: u32 }{
    .{ .name = "drug_like", .count = 0 },
    .{ .name = "adversarial", .count = 2000 },
};

/// Declares what the caller has arranged for running the other CPU
/// architecture's binaries. Hosts that need no help (same architecture, or
/// arm64 macOS running x86_64 under Rosetta) are detected instead of declared.
const ForeignExecutor = enum {
    /// Nothing arranged; the architecture axis fails unless the host can
    /// execute the foreign binaries on its own.
    none,
    /// qemu-user plus a cross libc, reached through `-fqemu --libc-runtimes`.
    qemu,
};

const OracleOptions = struct {
    target: std.Build.ResolvedTarget,
    optimize: std.builtin.OptimizeMode,
    oracle: *std.Build.Dependency,
    instrumented_fragment_builder: std.Build.LazyPath,
    instrumented_minimizer: std.Build.LazyPath,
};

const Oracle = struct {
    library: *std.Build.Step.Compile,
    abi: *std.Build.Step.Compile,
};

/// Builds the pinned upstream C++ for one target, plus the conformance-only
/// C ABI adapter over it. Both the native oracle and the other-architecture
/// oracle used by the stability classification come from here, so the two
/// differ in exactly one variable.
fn addOracle(b: *std.Build, options: OracleOptions) Oracle {
    const suffix = b.fmt("{s}-{s}", .{
        @tagName(options.target.result.cpu.arch),
        @tagName(options.target.result.os.tag),
    });

    // The upstream library is Boost-free and maeparser-free: only its CMake
    // test target needs them, which is why the harness and fixture reading
    // are re-hosted on the Zig side instead. Listing the pinned translation
    // units explicitly makes an upstream file appearing or disappearing a
    // build failure rather than a silent change in what the oracle contains.
    const library_module = b.createModule(.{
        .target = options.target,
        .optimize = options.optimize,
        .link_libcpp = true,
    });
    library_module.addIncludePath(b.path("conformance/include"));
    // The instrumented copies of two translation units are generated into the
    // cache, so their quoted includes no longer resolve next to the source.
    library_module.addIncludePath(options.oracle.path("."));
    library_module.addCSourceFiles(.{
        .root = options.oracle.path("."),
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
    library_module.addCSourceFile(.{
        .file = options.instrumented_fragment_builder,
        .flags = &.{ "-std=c++17", "-DIN_COORDGEN" },
    });
    library_module.addCSourceFile(.{
        .file = options.instrumented_minimizer,
        .flags = &.{ "-std=c++17", "-DIN_COORDGEN" },
    });
    library_module.addCSourceFile(.{
        .file = b.path("conformance/oracle_hook.cpp"),
        .flags = &.{ "-std=c++17", "-DIN_COORDGEN" },
    });
    const library = b.addLibrary(.{
        .name = b.fmt("coordgen-oracle-{s}", .{suffix}),
        .linkage = .static,
        .root_module = library_module,
    });

    // The adapter has the production-shaped C ABI plus the deliberately
    // broader probe API. It is conformance-only: neither artifact nor header
    // is installed with the native Zig library.
    const abi_module = b.createModule(.{
        .target = options.target,
        .optimize = options.optimize,
        .link_libcpp = true,
    });
    abi_module.addIncludePath(b.path("include"));
    abi_module.addIncludePath(b.path("conformance/include"));
    abi_module.addIncludePath(options.oracle.path("."));
    abi_module.addCSourceFile(.{
        .file = b.path("conformance/oracle_adapter.cpp"),
        .flags = &.{ "-std=c++17", "-Wall", "-Wextra", "-Werror" },
    });
    abi_module.linkLibrary(library);
    const abi = b.addLibrary(.{
        .name = b.fmt("coordgen-oracle-abi-{s}", .{suffix}),
        .linkage = .static,
        .root_module = abi_module,
    });

    return .{ .library = library, .abi = abi };
}

const CorpusRunnerOptions = struct {
    name: []const u8,
    target: std.Build.ResolvedTarget,
    optimize: std.builtin.OptimizeMode,
    conformance: *std.Build.Module,
    c_abi: *std.Build.Module,
    oracle_abi: *std.Build.Step.Compile,
    descending_allocator: bool,
};

fn addCorpusRunner(b: *std.Build, options: CorpusRunnerOptions) *std.Build.Step.Compile {
    const module = b.createModule(.{
        .root_source_file = b.path("tests/oracle_corpus_run.zig"),
        .target = options.target,
        .optimize = options.optimize,
        .imports = &.{
            .{ .name = "conformance", .module = options.conformance },
            .{ .name = "c_abi", .module = options.c_abi },
        },
        .link_libc = true,
        .link_libcpp = true,
    });
    if (options.descending_allocator) {
        module.addCSourceFile(.{
            .file = b.path("conformance/allocator_order.cpp"),
            .flags = &.{ "-std=c++17", "-Wall", "-Wextra", "-Werror" },
        });
    }
    module.linkLibrary(options.oracle_abi);
    return b.addExecutable(.{ .name = options.name, .root_module = module });
}

fn runCorpus(
    b: *std.Build,
    runner: *std.Build.Step.Compile,
    label: []const u8,
) [corpus_partitions.len]std.Build.LazyPath {
    var dumps: [corpus_partitions.len]std.Build.LazyPath = undefined;
    for (corpus_partitions, 0..) |partition, index| {
        const run = b.addRunArtifact(runner);
        run.expectExitCode(0);
        run.addArgs(&.{ "--partition", partition.name });
        run.addArgs(&.{ "--count", b.fmt("{d}", .{partition.count}) });
        run.addArg("--output");
        dumps[index] = run.addOutputFileArg(b.fmt("{s}-{s}.corpus", .{ label, partition.name }));
    }
    return dumps;
}

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

/// Every declared layer as named modules for one target, wired with exactly
/// the approved edges. A cross-compiled conformance harness needs its own
/// instance because a module belongs to one target.
const Layers = struct {
    core: *std.Build.Module,
    model: *std.Build.Module,
    geometry: *std.Build.Module,
    topology: *std.Build.Module,
    layout: *std.Build.Module,
    optimize_layer: *std.Build.Module,
    generator: *std.Build.Module,
    api: *std.Build.Module,
    c_abi: *std.Build.Module,
    conformance: *std.Build.Module,
};

fn createLayers(
    b: *std.Build,
    target: std.Build.ResolvedTarget,
    optimize: std.builtin.OptimizeMode,
) Layers {
    const layers: Layers = .{
        .core = createInternalModule(b, target, optimize, "src/core.zig"),
        .model = createInternalModule(b, target, optimize, "src/model.zig"),
        .geometry = createInternalModule(b, target, optimize, "src/geometry.zig"),
        .topology = createInternalModule(b, target, optimize, "build_support/empty_module.zig"),
        .layout = createInternalModule(b, target, optimize, "build_support/empty_module.zig"),
        .optimize_layer = createInternalModule(b, target, optimize, "build_support/empty_module.zig"),
        .generator = createInternalModule(b, target, optimize, "build_support/empty_module.zig"),
        .api = createInternalModule(b, target, optimize, "src/api.zig"),
        .c_abi = createInternalModule(b, target, optimize, "src/c_abi_types.zig"),
        .conformance = createInternalModule(b, target, optimize, "src/conformance.zig"),
    };
    const layered_modules = [_]*std.Build.Module{
        layers.core,
        layers.model,
        layers.geometry,
        layers.topology,
        layers.layout,
        layers.optimize_layer,
        layers.generator,
        layers.api,
        layers.c_abi,
        layers.conformance,
    };
    wireApprovedModuleEdges(&layered_modules);
    return layers;
}

pub fn build(b: *std.Build) !void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});

    const layers = createLayers(b, target, optimize);
    const core = layers.core;
    const model = layers.model;
    const geometry = layers.geometry;
    const api = layers.api;
    const c_abi = layers.c_abi;
    const conformance = layers.conformance;

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

    // Named-module imports mean a layer's tests belong to that layer's own
    // test binary: aggregating them behind a source-relative import in
    // src/coordgen.zig would silently stop running them.
    const layer_test_runs = blk: {
        const layer_modules = [_]struct { name: []const u8, module: *std.Build.Module }{
            .{ .name = "core-test", .module = core },
            .{ .name = "model-test", .module = model },
            .{ .name = "geometry-test", .module = geometry },
            .{ .name = "api-test", .module = api },
            .{ .name = "c-abi-test", .module = c_abi },
        };
        var runs: [layer_modules.len]*std.Build.Step.Run = undefined;
        for (layer_modules, 0..) |entry, index| {
            const tests = b.addTest(.{ .name = entry.name, .root_module = entry.module });
            runs[index] = b.addRunArtifact(tests);
        }
        break :blk runs;
    };
    const conformance_tests = b.addTest(.{
        .name = "conformance-module-test",
        .root_module = conformance,
    });
    const run_conformance_tests = b.addRunArtifact(conformance_tests);
    const corpus_classify_module = b.createModule(.{
        .root_source_file = b.path("tests/corpus_classify.zig"),
        .target = target,
        .optimize = optimize,
        .imports = &.{.{ .name = "conformance", .module = conformance }},
        .link_libc = false,
        .link_libcpp = false,
    });
    const corpus_classify_tests = b.addTest(.{
        .name = "corpus-classify-test",
        .root_module = corpus_classify_module,
    });
    const run_corpus_classify_tests = b.addRunArtifact(corpus_classify_tests);
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
    for (layer_test_runs) |run| test_step.dependOn(&run.step);
    test_step.dependOn(&run_conformance_tests.step);
    test_step.dependOn(&run_corpus_classify_tests.step);
    test_step.dependOn(&run_layer_tests.step);
    test_step.dependOn(&run_consumer_tests.step);

    const module_graph_step = b.step("module-graph-check", "Validate the approved named-module edge set");
    module_graph_step.dependOn(&run_layer_tests.step);
    // Zig's own compile errors catch an unapproved *named* import (the
    // module simply doesn't exist), but nothing stops a source-relative
    // import from reaching another layer's private files on disk, bypassing
    // the named-module table entirely. tools/check-module-imports is the
    // independent, CI-enforced check for that escape hatch; its --self-test
    // mode is the negative fixture proving it actually rejects one.
    const module_import_check = b.addSystemCommand(&.{ "python3", "tools/check-module-imports" });
    const module_import_self_test = b.addSystemCommand(&.{ "python3", "tools/check-module-imports", "--self-test" });
    module_graph_step.dependOn(&module_import_check.step);
    module_graph_step.dependOn(&module_import_self_test.step);

    const policy_command = b.addSystemCommand(&.{ "sh", "tools/check-build-policy" });
    const external_consumer = b.addRunFile(std.Build.LazyPath.zig_exe);
    external_consumer.addArgs(&.{ "build", "test", "--summary", "failures", "--cache-dir" });
    external_consumer.addDirectoryArg(std.Build.LazyPath.cache_root.path(b, "external-consumer"));
    external_consumer.setCwd(b.path("build_support/consumer"));
    const package_step = b.step("package-check", "Verify package and dependency isolation policy");
    package_step.dependOn(&policy_command.step);
    package_step.dependOn(&run_consumer_tests.step);
    package_step.dependOn(&external_consumer.step);
    // package-check is the gate every commit already runs; folding the
    // module-import escape-hatch check in here means it cannot be silently
    // skipped by only remembering `zig build test`.
    package_step.dependOn(module_graph_step);

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
    // The static asserts above check the header against itself; linking the
    // real library is what catches drift between the header's declared
    // coordgen_generate/coordgen_result_free signatures and what the Zig
    // implementation actually exports. A signature or name mismatch here is
    // a link error, not a silent pass.
    abi_layout_module.linkLibrary(library);
    const abi_layout = b.addExecutable(.{
        .name = "abi-layout-test",
        .root_module = abi_layout_module,
    });
    const run_abi_layout = b.addRunArtifact(abi_layout);
    run_abi_layout.expectExitCode(0);

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

    // `-fqemu` and `--libc-runtimes` are consumed by the build runner itself
    // and are not visible to this script: std.Build.enable_qemu is a vestigial
    // field that the 0.17-dev Maker never populates, so a build script that
    // reads it concludes "no executor" even when qemu is in use. The caller
    // therefore declares the executor separately. A false declaration is not a
    // way past the gate: the architecture-axis Run step still has to execute.
    const foreign_executor = b.option(
        ForeignExecutor,
        "foreign-executor",
        "How corpus-check runs other-architecture binaries the host cannot execute natively",
    ) orelse .none;

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

        const native_oracle = addOracle(b, .{
            .target = target,
            .optimize = optimize,
            .oracle = oracle,
            .instrumented_fragment_builder = instrumented_fragment_builder,
            .instrumented_minimizer = instrumented_minimizer,
        });
        const oracle_library = native_oracle.library;
        const oracle_abi = native_oracle.abi;
        // Deliberately never installed: the oracle exists for conformance
        // comparison only.
        oracle_step.dependOn(&oracle_library.step);
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

        // cgz-r07's sufficiency test: the RDKit-shaped flow, linked against the
        // pinned facade and run, rather than compiled for shape alone as
        // tests/abi_cpp_consumer.cpp is. The input-order contract is only
        // falsifiable when the facade is live, which is why this cannot live in
        // abi-check.
        const rdkit_consumer_module = b.createModule(.{
            .target = target,
            .optimize = optimize,
            .link_libcpp = true,
        });
        rdkit_consumer_module.addIncludePath(b.path("include"));
        rdkit_consumer_module.addCSourceFile(.{
            .file = b.path("tests/rdkit_consumer.cpp"),
            .flags = &.{ "-std=c++17", "-Wall", "-Wextra", "-Werror" },
        });
        rdkit_consumer_module.linkLibrary(oracle_abi);
        const rdkit_consumer = b.addExecutable(.{
            .name = "rdkit-consumer",
            .root_module = rdkit_consumer_module,
        });
        const run_rdkit_consumer = b.addRunArtifact(rdkit_consumer);
        run_rdkit_consumer.expectExitCode(0);
        oracle_step.dependOn(&run_rdkit_consumer.step);

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

        // Corpus stability classification. The same corpus runs through
        // oracle builds that differ in exactly one variable at a time:
        // architecture, or heap address order. Comparing the three dumps is
        // what assigns each input and observable its comparison tier.
        // The classification is a per-(architecture, toolchain,
        // optimize-mode) artifact, so the corpus oracles are pinned to one
        // optimize mode instead of following -Doptimize. That keeps a
        // classification comparable between runs, and keeps a 2000-member
        // corpus from being run through an unoptimized oracle three times.
        const corpus_optimize: std.builtin.OptimizeMode = .ReleaseFast;
        const corpus_layers = createLayers(b, target, corpus_optimize);
        const corpus_oracle = addOracle(b, .{
            .target = target,
            .optimize = corpus_optimize,
            .oracle = oracle,
            .instrumented_fragment_builder = instrumented_fragment_builder,
            .instrumented_minimizer = instrumented_minimizer,
        }).abi;
        const corpus_run = addCorpusRunner(b, .{
            .name = "oracle-corpus-run",
            .target = target,
            .optimize = corpus_optimize,
            .conformance = corpus_layers.conformance,
            .c_abi = corpus_layers.c_abi,
            .oracle_abi = corpus_oracle,
            .descending_allocator = false,
        });
        const corpus_run_descending = addCorpusRunner(b, .{
            .name = "oracle-corpus-run-descending",
            .target = target,
            .optimize = corpus_optimize,
            .conformance = corpus_layers.conformance,
            .c_abi = corpus_layers.c_abi,
            .oracle_abi = corpus_oracle,
            .descending_allocator = true,
        });
        const allocator_smoke_module = b.createModule(.{
            .target = target,
            .optimize = corpus_optimize,
            .link_libcpp = true,
        });
        allocator_smoke_module.addCSourceFiles(.{
            .files = &.{
                "conformance/allocator_order.cpp",
                "conformance/allocator_order_smoke.cpp",
            },
            .flags = &.{ "-std=c++17", "-Wall", "-Wextra", "-Werror" },
        });
        const allocator_smoke = b.addExecutable(.{
            .name = "allocator-order-smoke",
            .root_module = allocator_smoke_module,
        });
        const run_allocator_smoke = b.addRunArtifact(allocator_smoke);
        run_allocator_smoke.expectExitCode(0);

        // The architecture axis is a second oracle built for the other CPU
        // architecture of the same OS. Floating-point evaluation differences
        // there push the discrete DOF search into different local minima,
        // which is the divergence the classification has to measure rather
        // than assume.
        const cross_query: std.Target.Query = .{
            .cpu_arch = if (target.result.cpu.arch == .x86_64) .aarch64 else .x86_64,
            .os_tag = target.result.os.tag,
            .abi = target.result.abi,
        };
        const cross_target = b.resolveTargetQuery(cross_query);
        const cross_oracle = addOracle(b, .{
            .target = cross_target,
            .optimize = corpus_optimize,
            .oracle = oracle,
            .instrumented_fragment_builder = instrumented_fragment_builder,
            .instrumented_minimizer = instrumented_minimizer,
        }).abi;
        const cross_layers = createLayers(b, cross_target, corpus_optimize);
        const corpus_run_architecture = addCorpusRunner(b, .{
            .name = "oracle-corpus-run-architecture",
            .target = cross_target,
            .optimize = corpus_optimize,
            .conformance = cross_layers.conformance,
            .c_abi = cross_layers.c_abi,
            .oracle_abi = cross_oracle,
            .descending_allocator = false,
        });

        const baseline_dumps = runCorpus(b, corpus_run, "baseline");
        const descending_dumps = runCorpus(b, corpus_run_descending, "descending");
        const architecture_dumps = runCorpus(b, corpus_run_architecture, "architecture");

        // The drug-like partition is committed as tables in
        // src/conformance/corpus.zig; upstream's own SMILES parser is the
        // authority those tables are checked against.
        const smiles_reference_module = b.createModule(.{
            .target = target,
            .optimize = optimize,
            .link_libcpp = true,
        });
        smiles_reference_module.addIncludePath(oracle.path("."));
        smiles_reference_module.addCSourceFile(.{
            .file = b.path("conformance/smiles_reference.cpp"),
            .flags = &.{ "-std=c++17", "-Wall", "-Wextra", "-Werror" },
        });
        smiles_reference_module.linkLibrary(oracle_library);
        const smiles_reference = b.addExecutable(.{
            .name = "smiles-reference",
            .root_module = smiles_reference_module,
        });
        const run_smiles_reference = b.addRunArtifact(smiles_reference);
        const smiles_reference_dump = run_smiles_reference.captureStdOut(.{});

        const corpus_dump_module = b.createModule(.{
            .root_source_file = b.path("tests/corpus_dump.zig"),
            .target = target,
            .optimize = optimize,
            .imports = &.{.{ .name = "conformance", .module = conformance }},
            .link_libc = false,
            .link_libcpp = false,
        });
        const corpus_dump = b.addExecutable(.{
            .name = "corpus-dump",
            .root_module = corpus_dump_module,
        });
        const run_corpus_dump = b.addRunArtifact(corpus_dump);
        run_corpus_dump.addArgs(&.{ "--partition", "drug_like" });
        const corpus_dump_output = run_corpus_dump.captureStdOut(.{});

        const drug_like_diff = b.addSystemCommand(&.{ "diff", "-u" });
        drug_like_diff.addFileArg(smiles_reference_dump);
        drug_like_diff.addFileArg(corpus_dump_output);
        drug_like_diff.expectExitCode(0);

        const classify = b.addExecutable(.{
            .name = "corpus-classify",
            .root_module = corpus_classify_module,
        });
        const run_classify = b.addRunArtifact(classify);
        run_classify.expectExitCode(0);
        for (corpus_partitions, 0..) |partition, partition_index| {
            run_classify.addArgs(&.{ "--partition", partition.name });
            run_classify.addArg("--baseline");
            run_classify.addFileArg(baseline_dumps[partition_index]);
            run_classify.addArg("--architecture");
            run_classify.addFileArg(architecture_dumps[partition_index]);
            run_classify.addArg("--allocator-order");
            run_classify.addFileArg(descending_dumps[partition_index]);
        }
        run_classify.addArg("--expectations");
        run_classify.addFileArg(b.path("conformance/parity_expectations.tsv"));
        run_classify.addArg("--manifest");
        const generated_manifest = run_classify.addOutputFileArg("parity_manifest.tsv");

        // The manifest is published evidence about the build that produced
        // it, so refreshing it is an explicit act rather than a side effect
        // of running the gate.
        const publish_manifest = b.addUpdateSourceFiles();
        publish_manifest.addCopyFileToSource(generated_manifest, "conformance/parity_manifest.tsv");
        const manifest_step = b.step(
            "parity-manifest",
            "Regenerate the published parity manifest from this build",
        );
        manifest_step.dependOn(&publish_manifest.step);

        const corpus_step = b.step(
            "corpus-check",
            "Classify corpus stability against the recorded parity manifest",
        );
        corpus_step.dependOn(&drug_like_diff.step);
        corpus_step.dependOn(&run_classify.step);
        corpus_step.dependOn(&run_allocator_smoke.step);

        // A Run step whose binary the host cannot execute is skipped, and a
        // skipped step still leaves the build green. The architecture axis
        // would then quietly disappear from the classification, so its
        // absence is made a failure here instead.
        const host = b.graph.host.result;
        const cross = cross_target.result;
        const architecture_runnable = host.os.tag == cross.os.tag and
            (host.cpu.arch == cross.cpu.arch or
                // Rosetta 2 runs the x86_64 half on an arm64 macOS host.
                (host.os.tag == .macos and host.cpu.arch == .aarch64 and cross.cpu.arch == .x86_64) or
                (host.os.tag == .linux and foreign_executor == .qemu));
        if (!architecture_runnable) {
            const missing_executor = b.addFail(b.fmt(
                "corpus-check needs to run {s}-{s} binaries for the architecture axis; " ++
                    "on Linux install qemu-user and a matching cross libc, then pass " ++
                    "-Dforeign-executor=qemu -fqemu --libc-runtimes <prefix>; " ++
                    "otherwise run on a host that can execute them",
                .{ @tagName(cross.cpu.arch), @tagName(cross.os.tag) },
            ));
            corpus_step.dependOn(&missing_executor.step);
        }
        conformance_step.dependOn(corpus_step);

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
    const fuzz_pending = b.addFail("fuzz is not implemented yet; see cgz-7v2.4.4");
    fuzz_step.dependOn(&fuzz_pending.step);
}
