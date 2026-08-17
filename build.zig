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
    sanitize_c: std.zig.SanitizeC = .off,
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
        .sanitize_c = options.sanitize_c,
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
    // Repo-owned, so it is held to the strict flag set even though it is
    // compiled into the oracle library beside upstream's own sources.
    // -DIN_COORDGEN is what lets it see upstream internals; it is not a reason
    // to stop treating warnings as errors in code this repo wrote.
    library_module.addCSourceFile(.{
        .file = b.path("conformance/oracle_hook.cpp"),
        .flags = &.{ "-std=c++17", "-DIN_COORDGEN", "-Wall", "-Wextra", "-Werror" },
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
        .sanitize_c = options.sanitize_c,
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

/// One entry per declared layer, in `Layer` order. A layer is a *set* of
/// modules, not necessarily a single one: the first entry is the layer's
/// canonical module - the one other layers get when they import the layer by
/// name - and any further entry is a second module inside the same layer,
/// subject to the same row of `approved_edges` and additionally importing the
/// canonical module under its own layer name. The `c_abi` layer is the only
/// one that needs this today, and the reason is a linker constraint rather
/// than a size one; see src/c_abi/exports.zig.
fn wireApprovedModuleEdges(layer_modules: []const []const *std.Build.Module) void {
    layer_contract.validate() catch |err| {
        std.debug.panic("invalid approved module graph: {s}", .{@errorName(err)});
    };

    const layer_count = @typeInfo(layer_contract.Layer).@"enum".field_names.len;
    if (layer_modules.len != layer_count) @panic("build module table does not cover every declared layer");
    for (layer_modules) |members| {
        if (members.len == 0) @panic("a declared layer has no module");
    }

    // The build edge set is not a second hand-maintained list: every addImport
    // comes directly from the reviewed allow-list in src/module_layers.zig.
    // Every module of a layer gets that layer's outgoing edges, so splitting a
    // layer across files cannot quietly widen or narrow what it may reach.
    for (layer_contract.approved_edges) |edge| {
        const to = layer_modules[@backingInt(edge.to)][0];
        for (layer_modules[@backingInt(edge.from)]) |from| {
            from.addImport(@tagName(edge.to), to);
        }
    }

    // Intra-layer: a layer's secondary modules reach its canonical module by
    // the layer's own name. tools/check-module-imports permits that bare
    // import precisely because it does not cross a layer boundary.
    for (layer_modules, 0..) |members, index| {
        const layer: layer_contract.Layer = @fromBackingInt(@intCast(index));
        for (members[1..]) |member| member.addImport(@tagName(layer), members[0]);
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
    /// The c_abi layer's second module: the only Zig-side definition of the
    /// public coordgen_* symbols. Kept out of `c_abi` so that importing the
    /// ABI's types never drags a competing implementation into a link that
    /// already has one. Only the installed library's root should import it.
    c_abi_exports: *std.Build.Module,
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
        .topology = createInternalModule(b, target, optimize, "src/topology.zig"),
        .layout = createInternalModule(b, target, optimize, "src/layout.zig"),
        .optimize_layer = createInternalModule(b, target, optimize, "src/optimize/continuous.zig"),
        .generator = createInternalModule(b, target, optimize, "src/generator/minimal.zig"),
        .api = createInternalModule(b, target, optimize, "src/api.zig"),
        .c_abi = createInternalModule(b, target, optimize, "src/c_abi_types.zig"),
        .c_abi_exports = createInternalModule(b, target, optimize, "src/c_abi/exports.zig"),
        .conformance = createInternalModule(b, target, optimize, "src/conformance.zig"),
    };
    const layered_modules = [_][]const *std.Build.Module{
        &.{layers.core},
        &.{layers.model},
        &.{layers.geometry},
        &.{layers.topology},
        &.{layers.layout},
        &.{layers.optimize_layer},
        &.{layers.generator},
        &.{layers.api},
        &.{ layers.c_abi, layers.c_abi_exports },
        &.{layers.conformance},
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
    const topology = layers.topology;
    const generator = layers.generator;
    const api = layers.api;
    const c_abi = layers.c_abi;
    const c_abi_exports = layers.c_abi_exports;
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
            // The installed library is the one artifact that must contain the
            // public coordgen_* definitions, so it is the one artifact that
            // imports the exports module. Adding this import anywhere that
            // also links the oracle-backed ABI is a duplicate-symbol link
            // error, which is the point: cgz-r28.
            .{ .name = "c_abi_exports", .module = c_abi_exports },
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
    b.installFile("packaging/coordgen.pc", "lib/pkgconfig/coordgen.pc");

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
            .{ .name = "topology-test", .module = topology },
            .{ .name = "layout-test", .module = layers.layout },
            .{ .name = "optimize-test", .module = layers.optimize_layer },
            .{ .name = "generator-test", .module = generator },
            .{ .name = "api-test", .module = api },
            .{ .name = "c-abi-test", .module = c_abi },
            .{ .name = "c-abi-exports-test", .module = c_abi_exports },
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
    const native_minimal_module = b.createModule(.{
        .root_source_file = b.path("tests/native_minimal.zig"),
        .target = target,
        .optimize = optimize,
        .imports = &.{
            .{ .name = "api", .module = api },
            .{ .name = "generator", .module = generator },
        },
        .link_libc = false,
        .link_libcpp = false,
    });
    const native_minimal_tests = b.addTest(.{
        .name = "native-minimal-test",
        .root_module = native_minimal_module,
    });
    const run_native_minimal_tests = b.addRunArtifact(native_minimal_tests);
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
    test_step.dependOn(&run_native_minimal_tests.step);
    test_step.dependOn(&run_layer_tests.step);
    test_step.dependOn(&run_consumer_tests.step);

    // The parity ceiling is enumerated per corpus member, so it is checkable
    // without an oracle: regenerate each enumerated member and confirm it
    // still hashes to the bytes the row was measured against. Deliberately on
    // `test` rather than on `corpus-check`, because the point of cgz-r26 is
    // that the enumeration is portable and does not need the oracle build.
    const parity_ceiling_module = b.createModule(.{
        .root_source_file = b.path("tests/parity_ceiling_check.zig"),
        .target = target,
        .optimize = optimize,
        .imports = &.{.{ .name = "conformance", .module = conformance }},
        .link_libc = false,
        .link_libcpp = false,
    });
    const parity_ceiling_check = b.addExecutable(.{
        .name = "parity-ceiling-check",
        .root_module = parity_ceiling_module,
    });
    const run_parity_ceiling_check = b.addRunArtifact(parity_ceiling_check);
    run_parity_ceiling_check.expectExitCode(0);
    run_parity_ceiling_check.addArg("--ceiling");
    run_parity_ceiling_check.addFileArg(b.path("conformance/parity_ceiling.tsv"));
    const parity_ceiling_step = b.step(
        "parity-ceiling-check",
        "Verify the enumerated parity ceiling still names the members it was measured on",
    );
    parity_ceiling_step.dependOn(&run_parity_ceiling_check.step);
    test_step.dependOn(parity_ceiling_step);

    const coverage_step = b.step("coverage-check", "Validate requirements-to-test traceability");
    const coverage_check = b.addSystemCommand(&.{
        "python3",
        "tools/check-requirements-coverage.py",
        "--summary",
        "docs/architecture/REQUIREMENTS_COVERAGE.md",
    });
    const coverage_self_test = b.addSystemCommand(&.{
        "python3",
        "tools/check-requirements-coverage.py",
        "--self-test",
    });
    coverage_step.dependOn(&coverage_check.step);
    coverage_step.dependOn(&coverage_self_test.step);
    test_step.dependOn(coverage_step);

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
    // The install gate the cgz-r28 defect walked straight through, now in the
    // build graph rather than a documented two-command sequence. With no
    // PREFIX argument the tool installs into a scratch directory it creates
    // empty and destroys afterwards, so it can never certify files a previous
    // run left behind, and it compiles and runs examples/install_consumer.c
    // against that prefix alone - the check that would have caught an archive
    // exporting nothing.
    const install_isolation = b.addSystemCommand(&.{ "sh", "tools/check-install-isolation" });
    install_isolation.setEnvironmentVariable("CGZ_ZIG", b.graph.zig_exe);
    install_isolation.addPrefixedDirectoryArg(
        "--cache-dir=",
        std.Build.LazyPath.cache_root.path(b, "install-isolation"),
    );

    const package_step = b.step("package-check", "Verify package and dependency isolation policy");
    package_step.dependOn(&policy_command.step);
    package_step.dependOn(&run_consumer_tests.step);
    package_step.dependOn(&external_consumer.step);
    package_step.dependOn(&install_isolation.step);
    // The non-fuzzing backend list in tests/fuzz_backend.zig mirrors a private
    // declaration in the toolchain, so it can go stale in the direction that
    // fails open. This reads the switch out of the pinned toolchain in use and
    // requires the two to agree, continuously rather than only when someone
    // runs `zig build fuzz`.
    const fuzz_backend_mirror = b.addSystemCommand(&.{ "python3", "tools/check-fuzz-backend" });
    fuzz_backend_mirror.setEnvironmentVariable("CGZ_ZIG", b.graph.zig_exe);
    package_step.dependOn(&fuzz_backend_mirror.step);
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
    const sanitize_oracle = b.option(
        bool,
        "sanitize-oracle",
        "Instrument the native C++ oracle/facade with full C-family undefined-behavior checks",
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
    const regeneration_step = b.step("template-regeneration-check", "Verify pinned C++ and MAE template provenance");
    const performance_baseline_step = b.step(
        "performance-baseline",
        "Measure the absolute oracle baseline on the representative corpus",
    );
    const performance_step = b.step("performance-check", "Enforce native/oracle performance thresholds");
    const sanitizer_step = b.step("sanitizer-check", "Run oracle facade consumers with UB instrumentation");
    const performance_pending = b.addFail(
        "performance-check awaits the first native generation baseline and reviewed per-bucket ratios; see cgz-7v2.4.7",
    );
    performance_step.dependOn(&performance_pending.step);
    const template_fixture = b.createModule(.{
        .root_source_file = b.path("src/layout/templates/data.zig"),
        .target = target,
        .optimize = optimize,
    });
    const template_generator_module = b.createModule(.{
        .root_source_file = b.path("tests/template_generate.zig"),
        .target = target,
        .optimize = optimize,
        .imports = &.{
            .{ .name = "conformance", .module = conformance },
            .{ .name = "template_fixture", .module = template_fixture },
        },
    });
    const template_generator_tests = b.addTest(.{
        .name = "template-generate-test",
        .root_module = template_generator_module,
    });
    test_step.dependOn(&b.addRunArtifact(template_generator_tests).step);

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

        // Count public tests from their declarations, not from upstream's
        // CTest registration: its SMILES target accidentally runs
        // test_coordgen a second time. The inventory checker also pins the six
        // harness defects whose corrected local treatment is documented in
        // PARITY_MANIFEST.md.
        const public_test_inventory = b.addSystemCommand(&.{
            "python3",
            "tools/check-upstream-test-inventory.py",
        });
        public_test_inventory.addDirectoryArg(oracle.path("."));
        oracle_step.dependOn(&public_test_inventory.step);

        const upstream_coverage = b.addSystemCommand(&.{
            "python3",
            "tools/check-requirements-coverage.py",
            "--upstream",
        });
        upstream_coverage.addDirectoryArg(oracle.path("."));
        upstream_coverage.addArgs(&.{
            "--summary",
            "docs/architecture/REQUIREMENTS_COVERAGE.md",
        });
        oracle_step.dependOn(&upstream_coverage.step);

        const template_generator = b.addExecutable(.{
            .name = "template-generate",
            .root_module = template_generator_module,
        });

        const generated_runs = [_]*std.Build.Step.Run{
            b.addRunArtifact(template_generator),
            b.addRunArtifact(template_generator),
        };
        var generated: [2]std.Build.LazyPath = undefined;
        for (generated_runs, 0..) |run, index| {
            run.addFileArg(oracle.path("templates.mae"));
            generated[index] = run.addOutputFileArg(b.fmt("templates-normalized-{d}.zig", .{index}));
        }
        const compare_runs = b.addSystemCommand(&.{ "cmp", "-s" });
        compare_runs.addFileArg(generated[0]);
        compare_runs.addFileArg(generated[1]);
        const compare_fixture = b.addSystemCommand(&.{ "cmp", "-s" });
        compare_fixture.addFileArg(generated[0]);
        compare_fixture.addFileArg(b.path("src/layout/templates/data.zig"));

        const extract_cpp = b.addSystemCommand(&.{ "python3", "tools/extract-template-reference.py" });
        extract_cpp.addFileArg(oracle.path("CoordgenTemplates.cpp"));
        const extracted = extract_cpp.addOutputFileArg("templates-from-cpp.zig");
        const compare_cpp = b.addSystemCommand(&.{ "cmp", "-s" });
        compare_cpp.addFileArg(extracted);
        compare_cpp.addFileArg(b.path("src/layout/templates/data.zig"));

        const generate_raw = b.addRunArtifact(template_generator);
        generate_raw.addFileArg(oracle.path("templates.mae"));
        const mae_raw = generate_raw.addOutputFileArg("templates-mae-raw.zig");
        generate_raw.addArg("--raw");
        const extract_raw = b.addSystemCommand(&.{ "python3", "tools/extract-template-reference.py" });
        extract_raw.addFileArg(oracle.path("CoordgenTemplates.cpp"));
        const cpp_raw = extract_raw.addOutputFileArg("templates-cpp-raw.zig");
        extract_raw.addArg("--raw");
        const compare_raw = b.addSystemCommand(&.{ "cmp", "-s" });
        compare_raw.addFileArg(mae_raw);
        compare_raw.addFileArg(cpp_raw);
        regeneration_step.dependOn(&compare_runs.step);
        regeneration_step.dependOn(&compare_fixture.step);
        regeneration_step.dependOn(&compare_cpp.step);
        regeneration_step.dependOn(&compare_raw.step);

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
            .sanitize_c = if (sanitize_oracle) .full else .off,
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
            .sanitize_c = if (sanitize_oracle) .full else .off,
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
            .sanitize_c = if (sanitize_oracle) .full else .off,
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
            .sanitize_c = if (sanitize_oracle) .full else .off,
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
        if (sanitize_oracle) {
            sanitizer_step.dependOn(&run_oracle_abi_smoke.step);
            sanitizer_step.dependOn(&run_oracle_abi_cpp_smoke.step);
            sanitizer_step.dependOn(&run_rdkit_consumer.step);
        } else {
            const sanitizer_disabled = b.addFail(
                "sanitizer-check requires -Denable-oracle=true -Dsanitize-oracle=true",
            );
            sanitizer_step.dependOn(&sanitizer_disabled.step);
        }

        // Rehost all twenty public cases with repo-owned runners. The first
        // thirteen do not consume MAE fixtures. The fixture-backed seven read
        // deterministic dumps produced by our Zig MAE reader, so neither
        // runner depends on Boost.Test or maeparser.
        const public_tests_module = b.createModule(.{
            .target = target,
            .optimize = optimize,
            .link_libcpp = true,
        });
        public_tests_module.addIncludePath(oracle.path("."));
        public_tests_module.addCSourceFile(.{
            .file = b.path("conformance/public_test_rehost.cpp"),
            .flags = &.{ "-std=c++17", "-DIN_COORDGEN", "-Wall", "-Wextra", "-Werror" },
        });
        public_tests_module.linkLibrary(oracle_library);
        const public_tests = b.addExecutable(.{
            .name = "upstream-public-tests",
            .root_module = public_tests_module,
        });
        const run_public_tests = b.addRunArtifact(public_tests);
        run_public_tests.expectExitCode(0);
        oracle_step.dependOn(&run_public_tests.step);

        const mae_public_test_dump_module = b.createModule(.{
            .root_source_file = b.path("tests/mae_public_test_dump.zig"),
            .target = target,
            .optimize = optimize,
            .imports = &.{.{ .name = "conformance", .module = conformance }},
        });
        const mae_public_test_dump = b.addExecutable(.{
            .name = "mae-public-test-dump",
            .root_module = mae_public_test_dump_module,
        });
        const fixture_paths = [_][]const u8{
            "test/test.mae",
            "templates.mae",
            "test/testChirality.mae",
            "test/nonterminalMetalZobs.mae",
            "test/metalZobs.mae",
            "test/macrocycle.mae",
            "test/test_mol.mae",
        };
        var fixture_dumps: [fixture_paths.len]std.Build.LazyPath = undefined;
        for (fixture_paths, 0..) |fixture_path, index| {
            const dump_fixture = b.addRunArtifact(mae_public_test_dump);
            dump_fixture.addFileArg(oracle.path(fixture_path));
            fixture_dumps[index] = dump_fixture.addOutputFileArg(
                b.fmt("public-test-fixture-{d}.txt", .{index}),
            );
        }

        const fixture_public_tests_module = b.createModule(.{
            .target = target,
            .optimize = optimize,
            .link_libcpp = true,
        });
        fixture_public_tests_module.addIncludePath(oracle.path("."));
        fixture_public_tests_module.addCSourceFile(.{
            .file = b.path("conformance/public_fixture_test_rehost.cpp"),
            .flags = &.{ "-std=c++17", "-DIN_COORDGEN", "-Wall", "-Wextra", "-Werror" },
        });
        fixture_public_tests_module.linkLibrary(oracle_library);
        const fixture_public_tests = b.addExecutable(.{
            .name = "upstream-public-fixture-tests",
            .root_module = fixture_public_tests_module,
        });
        const run_fixture_public_tests = b.addRunArtifact(fixture_public_tests);
        for (fixture_dumps) |dump| {
            run_fixture_public_tests.addFileArg(dump);
        }
        run_fixture_public_tests.expectExitCode(0);
        oracle_step.dependOn(&run_fixture_public_tests.step);

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
        const benchmark_module = b.createModule(.{
            .root_source_file = b.path("tests/oracle_benchmark.zig"),
            .target = target,
            .optimize = corpus_optimize,
            .imports = &.{
                .{ .name = "conformance", .module = corpus_layers.conformance },
                .{ .name = "c_abi", .module = corpus_layers.c_abi },
            },
            .link_libc = true,
            .link_libcpp = true,
        });
        benchmark_module.addCSourceFile(.{
            .file = b.path("conformance/benchmark_clock.c"),
            .flags = &.{ "-std=c11", "-Wall", "-Wextra", "-Werror", "-pedantic" },
        });
        benchmark_module.linkLibrary(corpus_oracle);
        const benchmark = b.addExecutable(.{
            .name = "oracle-performance-baseline",
            .root_module = benchmark_module,
        });
        const benchmark_tests = b.addTest(.{
            .name = "oracle-performance-contract-test",
            .root_module = benchmark_module,
        });
        oracle_step.dependOn(&b.addRunArtifact(benchmark_tests).step);
        const run_benchmark = b.addRunArtifact(benchmark);
        // Timings are evidence, not deterministic cache outputs: print them
        // for capture by the caller and rerun the process on every invocation.
        run_benchmark.stdio = .inherit;
        performance_baseline_step.dependOn(&run_benchmark.step);
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
        // The fraction ceilings bound how many pairs may be order-unstable;
        // this enumerates which. Both are checked, and an order-unstable pair
        // missing from the enumeration fails closed (cgz-r13, cgz-r26).
        run_classify.addArg("--ceiling");
        run_classify.addFileArg(b.path("conformance/parity_ceiling.tsv"));
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

        const native_macrocycle_module = b.createModule(.{
            .root_source_file = b.path("tests/native_macrocycle_mae.zig"),
            .target = target,
            .optimize = optimize,
            .imports = &.{
                .{ .name = "api", .module = api },
                .{ .name = "conformance", .module = conformance },
                .{ .name = "generator", .module = generator },
            },
            .link_libc = false,
            .link_libcpp = false,
        });
        const native_macrocycle = b.addExecutable(.{
            .name = "native-macrocycle-mae",
            .root_module = native_macrocycle_module,
        });
        const run_native_macrocycle = b.addRunArtifact(native_macrocycle);
        run_native_macrocycle.expectExitCode(0);
        run_native_macrocycle.addFileArg(oracle.path("test/macrocycle.mae"));
        oracle_step.dependOn(&run_native_macrocycle.step);

        // The differential runner (cgz-7v2.4.2). It imports `api` for the
        // native side and links the oracle's C ABI for the other, which is
        // only sound because the native C exports are absent from this module
        // - see tests/native_oracle_diff.zig and cgz-r28.
        const native_diff_module = b.createModule(.{
            .root_source_file = b.path("tests/native_oracle_diff.zig"),
            .target = target,
            .optimize = optimize,
            .imports = &.{
                .{ .name = "api", .module = api },
                .{ .name = "core", .module = core },
                .{ .name = "conformance", .module = conformance },
                .{ .name = "c_abi", .module = c_abi },
            },
            .link_libc = true,
            .link_libcpp = true,
        });
        native_diff_module.linkLibrary(oracle_abi);
        const native_diff = b.addExecutable(.{
            .name = "native-oracle-diff",
            .root_module = native_diff_module,
        });
        const run_native_diff = b.addRunArtifact(native_diff);
        run_native_diff.expectExitCode(0);
        // drug_like first: it is the partition the Phase 3 gate in
        // SUCCESS_CRITERIA.md is stated over, and every member is
        // coordinate-stable on both axes, so a mismatch here is the port's.
        run_native_diff.addArgs(&.{ "--partition", "drug_like" });
        run_native_diff.addArg("--ceiling");
        run_native_diff.addFileArg(b.path("conformance/parity_ceiling.tsv"));
        run_native_diff.addArg("--baseline");
        const native_diff_baseline = run_native_diff.addOutputFileArg("native_baseline.tsv");
        const publish_baseline = b.addUpdateSourceFiles();
        publish_baseline.addCopyFileToSource(native_diff_baseline, "conformance/native_baseline.tsv");
        const native_diff_step = b.step(
            "native-diff",
            "Compare native against the pinned oracle over a corpus partition",
        );
        native_diff_step.dependOn(&run_native_diff.step);
        const native_baseline_step = b.step(
            "native-baseline",
            "Regenerate the published native-vs-oracle baseline from this build",
        );
        native_baseline_step.dependOn(&publish_baseline.step);
        conformance_step.dependOn(native_diff_step);

        conformance_step.dependOn(oracle_step);
    } else {
        const oracle_disabled = b.addFail(
            "oracle steps require -Denable-oracle=true; ordinary builds never fetch conformance sources",
        );
        oracle_step.dependOn(&oracle_disabled.step);
        conformance_step.dependOn(&oracle_disabled.step);
        regeneration_step.dependOn(&oracle_disabled.step);
        performance_baseline_step.dependOn(&oracle_disabled.step);
        sanitizer_step.dependOn(&oracle_disabled.step);
    }

    // Both examples execute as real consumers: Zig reads its source from the
    // external dependency package, while C compiles against a fresh scratch
    // install prefix. Sharing these runs with package-check keeps the package
    // and examples gates from drifting apart.
    const examples_step = b.step("examples", "Run package and installed-artifact examples");
    examples_step.dependOn(&external_consumer.step);
    examples_step.dependOn(&install_isolation.step);

    const fuzz_step = b.step("fuzz", "Run the platform-selected fuzz harness when available");
    const fuzz_pending = b.addFail("fuzz is not implemented yet; see cgz-7v2.4.4");
    fuzz_step.dependOn(&fuzz_pending.step);
    // Attached now, before the harness exists, so cgz-7v2.4.4 inherits it
    // rather than having to remember it: on a compiler backend whose
    // std.testing.fuzz returns immediately, every fuzz target is a no-op that
    // reports success. This makes `zig build fuzz` fail by name instead. It is
    // deliberately not on `test` - a backend that cannot fuzz is no reason to
    // block unrelated work. See cgz-r27.
    const fuzz_backend_module = b.createModule(.{
        .root_source_file = b.path("tests/fuzz_backend.zig"),
        .target = target,
        .optimize = optimize,
    });
    const fuzz_backend_tests = b.addTest(.{
        .name = "fuzz-backend-test",
        .root_module = fuzz_backend_module,
    });
    fuzz_step.dependOn(&b.addRunArtifact(fuzz_backend_tests).step);

    // A step with no dependencies succeeds instantly and reports success, which
    // is the cgz-r20 failure mode. tools/check-gate-strength used to assert
    // this by grepping build.zig for the string "<step>.dependOn(" - a text
    // property, satisfiable by a call inside a branch that never runs. This
    // asks the constructed graph instead, which is the thing that has to be
    // true, and it runs at configure time on every single invocation rather
    // than only under package-check.
    assertNoFalseGreenSteps(&.{
        test_step,
        module_graph_step,
        package_step,
        abi_step,
        oracle_step,
        conformance_step,
        coverage_step,
        examples_step,
        regeneration_step,
        performance_baseline_step,
        performance_step,
        sanitizer_step,
        fuzz_step,
    });
}

/// Every publicly reachable step must do something. Configure-time panic
/// rather than a build error because a gate that certifies nothing is a
/// defect in the build description itself, not in what is being built.
fn assertNoFalseGreenSteps(steps: []const *std.Build.Step) void {
    for (steps) |step| {
        if (step.dependencies.items.len == 0) {
            std.debug.panic(
                "build step '{s}' has no dependencies and would report success without checking anything",
                .{step.name},
            );
        }
    }
}
