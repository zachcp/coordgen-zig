//! Typed comparison semantics for the future native-vs-oracle runner.
//!
//! A comparison mode is selected for every input *and* observable. Selecting
//! one mode for a whole molecule would hide stable structural evidence merely
//! because that molecule's coordinates are unstable.

const std = @import("std");

pub const Tier = enum { exact, tolerant, invariant_statistical };

pub const Observable = enum {
    status,
    clean_pose,
    coordinates,
    input_to_internal,
    internal_to_input,
    effective_bond_orders,
    bond_displays,
    atom_stereo,
    probe_status,
    probe_clean_pose,
    morgan_ranks,
    rings,
    rings_set,
    fragments,
    fragments_set,
    dofs,
    dofs_set,
    dof_penalties,
    dof_penalties_set,
    template_mappings,
    template_mappings_set,
    components,
    components_set,
    component_transforms,
    component_transforms_set,
};

pub const Stability = struct {
    architecture: bool,
    allocator_order: bool,

    pub fn exact(self: Stability) bool {
        return self.architecture and self.allocator_order;
    }
};

/// The tier an observable is entitled to when nothing perturbs it: structural
/// values are exact, coordinates and derived floats are tolerant. Stability is
/// not an input here, so this is the ceiling a pair can never exceed.
pub fn baseTier(observable: Observable) Tier {
    return switch (observable) {
        .coordinates,
        .dof_penalties,
        .dof_penalties_set,
        .component_transforms,
        .component_transforms_set,
        => .tolerant,
        else => .exact,
    };
}

/// Structural values are exact whenever their own stability record is exact.
/// Coordinates and derived floats use normalized/tolerant comparisons only if
/// the oracle is stable, otherwise they use the aggregate quality tier.
///
/// **Superseded for the differential runner; use `differentialComparison`.**
/// This reads `Stability.exact()`, which is both axes, so it predates
/// `cgz-r13`'s finding that the runner never observes the architecture axis
/// (both sides share a build). Feeding it a member that is arch-unstable but
/// order-stable yields `.invariant_statistical`, and `cgz-r30` records that
/// the invariant tier is unpopulated at the pin — so that answer is wrong for
/// any real comparison.
///
/// Kept rather than deleted because `cgz-r06` froze this tier vocabulary and
/// the test below encodes its meaning; nothing on a live path calls it.
pub fn tierFor(observable: Observable, stability: Stability) Tier {
    if (!stability.exact()) return .invariant_statistical;
    return baseTier(observable);
}

pub const bond_length: f64 = 50;

// ---------------------------------------------------------------------------
// The portable parity-ceiling enumeration (cgz-r26)
// ---------------------------------------------------------------------------

/// Deviation at which a coordinate difference stops being float jitter and
/// becomes a different layout, taken from the measured bimodal distribution in
/// SUCCESS_CRITERIA.md: only 2.3% of divergent members land in the [0.1, 1.0)
/// valley, and 23.6% sit at or above 1.0. A per-pair ceiling at or above this
/// boundary would excuse a relayout, which `cgz-r13` forbids, so the loader
/// rejects one.
pub const relayout_boundary_bond_lengths: f64 = 1.0;

pub const CeilingParseError = error{
    MalformedCeilingRow,
    UnknownCeilingObservable,
    ExactObservableAtCeiling,
    CeilingBoundOnNonCoordinateObservable,
    CeilingBoundNotPositive,
    CeilingBoundReachesRelayoutBoundary,
    MalformedCeilingDigest,
    MissingCeilingEvidence,
    DuplicateCeilingRow,
};

/// One enumerated (member, observable) pair at the pointer-order parity
/// ceiling. `member` is the portable corpus identity `{partition}/{index}`;
/// `input_sha256` pins the bytes that identity generates so the identity keeps
/// meaning the same molecule across a corpus regeneration.
pub const CeilingRow = struct {
    member: []const u8,
    observable: Observable,
    /// `null` means the pair keeps its normal tier's tolerance: being at the
    /// ceiling excuses the pair from *exact* comparison, it does not widen
    /// anything. A value replaces that tolerance for this pair alone.
    order_ceiling_bond_lengths: ?f64,
    input_sha256: []const u8,
    evidence: []const u8,
};

/// A borrowed view over `conformance/parity_ceiling.tsv`. The file is small
/// and read once, so lookup is a linear scan and the table allocates nothing.
pub const CeilingTable = struct {
    text: []const u8,

    pub fn iterate(self: CeilingTable) RowIterator {
        return .{ .lines = std.mem.splitScalar(u8, self.text, '\n') };
    }

    /// The row governing this pair, or null when the pair is not enumerated.
    /// A null answer is a failure at every call site: an unenumerated pair
    /// that measured order-unstable is exactly what the enumeration exists to
    /// reject.
    pub fn find(
        self: CeilingTable,
        member: []const u8,
        observable: Observable,
    ) CeilingParseError!?CeilingRow {
        var rows = self.iterate();
        while (try rows.next()) |row| {
            if (row.observable == observable and std.mem.eql(u8, row.member, member)) {
                return row;
            }
        }
        return null;
    }

    /// Rejects a file that repeats a pair, which would otherwise let one row
    /// shadow a stricter one depending on scan order.
    pub fn validateUnique(self: CeilingTable) CeilingParseError!void {
        var outer = self.iterate();
        var position: usize = 0;
        while (try outer.next()) |row| : (position += 1) {
            var inner = self.iterate();
            var other_position: usize = 0;
            while (try inner.next()) |other| : (other_position += 1) {
                if (other_position >= position) break;
                if (other.observable == row.observable and
                    std.mem.eql(u8, other.member, row.member))
                {
                    return error.DuplicateCeilingRow;
                }
            }
        }
    }

    pub const RowIterator = struct {
        lines: std.mem.SplitIterator(u8, .scalar),

        pub fn next(self: *RowIterator) CeilingParseError!?CeilingRow {
            while (self.lines.next()) |line| {
                const trimmed = std.mem.trim(u8, line, " \t\r");
                if (trimmed.len == 0 or trimmed[0] == '#') continue;
                return try parseCeilingRow(trimmed);
            }
            return null;
        }
    };
};

fn parseCeilingRow(line: []const u8) CeilingParseError!CeilingRow {
    var fields = std.mem.splitScalar(u8, line, '\t');
    const member = fields.next() orelse return error.MalformedCeilingRow;
    const observable_name = fields.next() orelse return error.MalformedCeilingRow;
    const bound_text = fields.next() orelse return error.MalformedCeilingRow;
    const digest = fields.next() orelse return error.MalformedCeilingRow;
    const evidence = fields.next() orelse return error.MalformedCeilingRow;
    if (fields.next() != null) return error.MalformedCeilingRow;
    if (member.len == 0 or std.mem.indexOfScalar(u8, member, '/') == null) {
        return error.MalformedCeilingRow;
    }
    if (evidence.len == 0) return error.MissingCeilingEvidence;

    const observable = std.meta.stringToEnum(Observable, observable_name) orelse
        return error.UnknownCeilingObservable;
    // cgz-r13: structural observables stay exact requirements on ceiling
    // members too. Admitting one here would be that weakening, so it is not
    // representable rather than merely discouraged.
    if (baseTier(observable) == .exact) return error.ExactObservableAtCeiling;

    var bound: ?f64 = null;
    if (!std.mem.eql(u8, bound_text, "-")) {
        // Only `coordinates` carries a published per-member magnitude; the
        // classifier measures no deviation for any other observable, so a
        // number on one would be invented rather than measured.
        if (observable != .coordinates) return error.CeilingBoundOnNonCoordinateObservable;
        const value = std.fmt.parseFloat(f64, bound_text) catch
            return error.MalformedCeilingRow;
        if (!std.math.isFinite(value) or value <= 0) return error.CeilingBoundNotPositive;
        if (value >= relayout_boundary_bond_lengths) {
            return error.CeilingBoundReachesRelayoutBoundary;
        }
        bound = value;
    }

    if (digest.len != 64) return error.MalformedCeilingDigest;
    for (digest) |character| {
        const hexadecimal = (character >= '0' and character <= '9') or
            (character >= 'a' and character <= 'f');
        if (!hexadecimal) return error.MalformedCeilingDigest;
    }

    return .{
        .member = member,
        .observable = observable,
        .order_ceiling_bond_lengths = bound,
        .input_sha256 = digest,
        .evidence = evidence,
    };
}

pub const Comparison = struct {
    tier: Tier,
    /// Zero for the exact tier. For the tolerant tier this is the tolerance the
    /// runner must use for this pair — the global `T` unless the enumeration
    /// published a wider per-pair bound.
    tolerance_bond_lengths: f64,
};

pub const DifferentialError = error{UnenumeratedOrderInstability};

/// The differential runner's per-pair rule.
///
/// `order_stable` is the only stability input, and deliberately so: the runner
/// executes oracle and native in one process, build, and target
/// (`RunIdentity.requireSameBuild`), so the architecture axis is not a
/// difference it can observe. `cgz-r13` scopes the ceiling to the heap-address
/// order axis for that reason.
///
/// A pair that is order-unstable with no enumerated row fails. Falling through
/// to the invariant tier instead would let a real defect hide behind an axis
/// nobody enumerated, which is the budget framing `cgz-r13` rejected.
pub fn differentialComparison(
    observable: Observable,
    order_stable: bool,
    ceiling: ?CeilingRow,
    tolerance_bond_lengths: f64,
) DifferentialError!Comparison {
    const base = baseTier(observable);
    if (order_stable) return .{
        .tier = base,
        .tolerance_bond_lengths = if (base == .exact) 0 else tolerance_bond_lengths,
    };
    const row = ceiling orelse return error.UnenumeratedOrderInstability;
    if (row.observable != observable) return error.UnenumeratedOrderInstability;
    return .{
        .tier = .tolerant,
        .tolerance_bond_lengths = row.order_ceiling_bond_lengths orelse tolerance_bond_lengths,
    };
}

pub const Point = struct { x: f64, y: f64 };

/// Reflection is unavailable unless the fixture records an achirality proof.
/// The default forbids it, preventing a mirrored chiral layout from passing a
/// translation/rotation-normalized comparison.
pub const ReflectionPolicy = struct {
    achiral_proof: ?[]const u8 = null,

    pub fn allowsReflection(self: ReflectionPolicy) bool {
        return if (self.achiral_proof) |proof| proof.len != 0 else false;
    }
};

pub const CoordinateResult = struct {
    matches: bool,
    max_deviation_bond_lengths: f64,
};

/// Compare indexed 2D coordinates after translation and rotation. A second,
/// reflected alignment is attempted only when the fixture opts in with proof.
pub fn compareCoordinates(
    reference: []const Point,
    candidate: []const Point,
    tolerance_bond_lengths: f64,
    reflection: ReflectionPolicy,
) !CoordinateResult {
    if (reference.len != candidate.len) return error.CountMismatch;
    if (!std.math.isFinite(tolerance_bond_lengths) or tolerance_bond_lengths < 0) {
        return error.InvalidTolerance;
    }
    const direct = try alignedDeviation(reference, candidate, false);
    if (direct <= tolerance_bond_lengths) return .{
        .matches = true,
        .max_deviation_bond_lengths = direct,
    };
    if (!reflection.allowsReflection()) return .{
        .matches = false,
        .max_deviation_bond_lengths = direct,
    };
    const reflected_deviation = try alignedDeviation(reference, candidate, true);
    return .{
        .matches = reflected_deviation <= tolerance_bond_lengths,
        .max_deviation_bond_lengths = reflected_deviation,
    };
}

fn alignedDeviation(reference: []const Point, candidate: []const Point, reflect: bool) !f64 {
    if (reference.len == 0) return 0;
    if (reference.len == 1) return distance(reference[0], candidate[0]) / bond_length;

    const candidate_first = reflected(candidate[0], reflect);
    var anchor: ?usize = null;
    for (reference[1..], candidate[1..], 1..) |reference_point, candidate_point, index| {
        if (norm(subtract(reference_point, reference[0])) != 0 and
            norm(subtract(reflected(candidate_point, reflect), candidate_first)) != 0)
        {
            anchor = index;
            break;
        }
    }
    const anchor_index = anchor orelse return error.DegenerateAlignment;
    const reference_vector = subtract(reference[anchor_index], reference[0]);
    const candidate_vector = subtract(reflected(candidate[anchor_index], reflect), candidate_first);
    const denominator = norm(reference_vector) * norm(candidate_vector);
    const cosine = dot(candidate_vector, reference_vector) / denominator;
    const sine = cross(candidate_vector, reference_vector) / denominator;

    var max_deviation: f64 = 0;
    for (reference, candidate) |reference_point, candidate_point| {
        const translated = subtract(reflected(candidate_point, reflect), candidate_first);
        const rotated = Point{
            .x = translated.x * cosine - translated.y * sine + reference[0].x,
            .y = translated.x * sine + translated.y * cosine + reference[0].y,
        };
        max_deviation = @max(max_deviation, distance(reference_point, rotated) / bond_length);
    }
    return max_deviation;
}

fn reflected(point: Point, enabled: bool) Point {
    return if (enabled) .{ .x = point.x, .y = -point.y } else point;
}

fn subtract(left: Point, right: Point) Point {
    return .{ .x = left.x - right.x, .y = left.y - right.y };
}

fn norm(point: Point) f64 {
    return std.math.sqrt(dot(point, point));
}

fn dot(left: Point, right: Point) f64 {
    return left.x * right.x + left.y * right.y;
}

fn cross(left: Point, right: Point) f64 {
    return left.x * right.y - left.y * right.x;
}

fn distance(left: Point, right: Point) f64 {
    return norm(subtract(left, right));
}

/// Coordinate parity is valid only within one process/build/target. This
/// identity belongs to each baseline instead of committing cross-platform
/// coordinate goldens that would create false failures.
pub const RunIdentity = struct {
    process_id: []const u8,
    target: []const u8,
    toolchain: []const u8,
    optimize_mode: []const u8,

    pub fn requireSameBuild(oracle: RunIdentity, native: RunIdentity) !void {
        if (!std.mem.eql(u8, oracle.process_id, native.process_id) or
            !std.mem.eql(u8, oracle.target, native.target) or
            !std.mem.eql(u8, oracle.toolchain, native.toolchain) or
            !std.mem.eql(u8, oracle.optimize_mode, native.optimize_mode))
        {
            return error.MismatchedRunIdentity;
        }
    }
};

/// Lower values are better. The invariant tier aggregates these per corpus
/// partition and requires native to be no worse than oracle plus a named
/// margin for every metric.
pub const QualityMetrics = struct {
    clash_score: f64,
    bond_length_rms: f64,
    bond_angle_deviation: f64,
    bond_crossings: u32,
    atoms_inside_rings: u32,
};

pub const QualityMargins = struct {
    clash_score: f64 = 0,
    bond_length_rms: f64 = 0,
    bond_angle_deviation: f64 = 0,
    bond_crossings: u32 = 0,
    atoms_inside_rings: u32 = 0,
};

pub fn noWorseThan(oracle: QualityMetrics, native: QualityMetrics, margin: QualityMargins) bool {
    return std.math.isFinite(oracle.clash_score) and
        std.math.isFinite(oracle.bond_length_rms) and
        std.math.isFinite(oracle.bond_angle_deviation) and
        std.math.isFinite(native.clash_score) and
        std.math.isFinite(native.bond_length_rms) and
        std.math.isFinite(native.bond_angle_deviation) and
        std.math.isFinite(margin.clash_score) and
        std.math.isFinite(margin.bond_length_rms) and
        std.math.isFinite(margin.bond_angle_deviation) and
        margin.clash_score >= 0 and
        margin.bond_length_rms >= 0 and
        margin.bond_angle_deviation >= 0 and
        native.clash_score <= oracle.clash_score + margin.clash_score and
        native.bond_length_rms <= oracle.bond_length_rms + margin.bond_length_rms and
        native.bond_angle_deviation <= oracle.bond_angle_deviation + margin.bond_angle_deviation and
        native.bond_crossings <= oracle.bond_crossings +| margin.bond_crossings and
        native.atoms_inside_rings <= oracle.atoms_inside_rings +| margin.atoms_inside_rings;
}

test "tier selection stays per input and observable" {
    const stable = Stability{ .architecture = true, .allocator_order = true };
    const unstable_coordinates = Stability{ .architecture = false, .allocator_order = true };
    try std.testing.expectEqual(Tier.exact, tierFor(.rings, stable));
    try std.testing.expectEqual(Tier.tolerant, tierFor(.coordinates, stable));
    try std.testing.expectEqual(Tier.invariant_statistical, tierFor(.coordinates, unstable_coordinates));
}

/// Digest-shaped filler: these fixtures exercise parsing and tier
/// selection, never identity, which tests/parity_ceiling_check.zig owns.
const test_digest = "0000000000000000000000000000000000000000000000000000000000000000";

const ceiling_fixture =
    "# member\tobservable\torder_ceiling_bl\tinput_sha256\tevidence\n" ++
    "adversarial/917\tcomponent_transforms\t-\t" ++ test_digest ++ "\tcgz-r13\n" ++
    "adversarial/1695\tcoordinates\t0.75\t" ++ test_digest ++ "\tcgz-r13\n";

test "the ceiling table parses, looks up, and rejects a repeated pair" {
    const table: CeilingTable = .{ .text = ceiling_fixture };
    try table.validateUnique();

    const found = (try table.find("adversarial/1695", .coordinates)).?;
    try std.testing.expectEqual(@as(?f64, 0.75), found.order_ceiling_bond_lengths);
    try std.testing.expectEqualStrings("cgz-r13", found.evidence);

    const unwidened = (try table.find("adversarial/917", .component_transforms)).?;
    try std.testing.expectEqual(@as(?f64, null), unwidened.order_ceiling_bond_lengths);

    try std.testing.expectEqual(
        @as(?CeilingRow, null),
        try table.find("adversarial/1695", .component_transforms),
    );

    const repeated: CeilingTable = .{ .text = ceiling_fixture ++ ceiling_fixture };
    try std.testing.expectError(error.DuplicateCeilingRow, repeated.validateUnique());
}

test "the ceiling file cannot express a weakening" {
    // A structural observable may never be admitted: cgz-r13 keeps those exact
    // on ceiling members too.
    try std.testing.expectError(error.ExactObservableAtCeiling, parseCeilingRow(
        "adversarial/1695\trings\t-\t" ++ test_digest ++ "\tcgz-r13",
    ));
    // A bound at or above the relayout boundary would excuse a different
    // layout rather than float jitter.
    try std.testing.expectError(error.CeilingBoundReachesRelayoutBoundary, parseCeilingRow(
        "adversarial/1695\tcoordinates\t1.0\t" ++ test_digest ++ "\tcgz-r13",
    ));
    // No observable but `coordinates` has a published magnitude to bound.
    try std.testing.expectError(error.CeilingBoundOnNonCoordinateObservable, parseCeilingRow(
        "adversarial/917\tcomponent_transforms\t0.5\t" ++ test_digest ++ "\tcgz-r13",
    ));
    // Identity, provenance, and pinned input bytes are all mandatory.
    try std.testing.expectError(error.MalformedCeilingDigest, parseCeilingRow(
        "adversarial/1695\tcoordinates\t0.75\tshort\tcgz-r13",
    ));
    try std.testing.expectError(error.MissingCeilingEvidence, parseCeilingRow(
        "adversarial/1695\tcoordinates\t0.75\t" ++ test_digest ++ "\t",
    ));
    try std.testing.expectError(error.MalformedCeilingRow, parseCeilingRow(
        "adversarial/1695\tcoordinates\t0.75\t" ++ test_digest,
    ));
}

test "an order-unstable pair with no enumerated row fails closed" {
    const table: CeilingTable = .{ .text = ceiling_fixture };
    const tolerance: f64 = 0.1;

    // Stable pairs keep their ordinary tier and the global tolerance.
    const stable_rings = try differentialComparison(.rings, true, null, tolerance);
    try std.testing.expectEqual(Tier.exact, stable_rings.tier);
    try std.testing.expectEqual(@as(f64, 0), stable_rings.tolerance_bond_lengths);

    // Enumerated with a published bound: tolerant against that bound alone.
    const widened = try differentialComparison(
        .coordinates,
        false,
        try table.find("adversarial/1695", .coordinates),
        tolerance,
    );
    try std.testing.expectEqual(Tier.tolerant, widened.tier);
    try std.testing.expectEqual(@as(f64, 0.75), widened.tolerance_bond_lengths);

    // Enumerated without one: excused from exact comparison, nothing widened.
    const unwidened = try differentialComparison(
        .component_transforms,
        false,
        try table.find("adversarial/917", .component_transforms),
        tolerance,
    );
    try std.testing.expectEqual(Tier.tolerant, unwidened.tier);
    try std.testing.expectEqual(tolerance, unwidened.tolerance_bond_lengths);

    // Not enumerated, and a row for a different observable does not stand in.
    try std.testing.expectError(
        error.UnenumeratedOrderInstability,
        differentialComparison(.coordinates, false, null, tolerance),
    );
    try std.testing.expectError(
        error.UnenumeratedOrderInstability,
        differentialComparison(
            .component_transforms,
            false,
            try table.find("adversarial/1695", .coordinates),
            tolerance,
        ),
    );
}

test "translation and rotation pass while reflection fails by default" {
    const reference = [_]Point{ .{ .x = 0, .y = 0 }, .{ .x = 50, .y = 0 }, .{ .x = 12.5, .y = 37.5 } };
    const rotated = [_]Point{ .{ .x = 200, .y = 300 }, .{ .x = 200, .y = 350 }, .{ .x = 162.5, .y = 312.5 } };
    const mirrored = [_]Point{ .{ .x = 0, .y = 0 }, .{ .x = 50, .y = 0 }, .{ .x = 12.5, .y = -37.5 } };

    try std.testing.expect((try compareCoordinates(&reference, &rotated, 1e-12, .{})).matches);
    try std.testing.expect(!(try compareCoordinates(&reference, &mirrored, 1e-12, .{})).matches);
    try std.testing.expect((try compareCoordinates(&reference, &mirrored, 1e-12, .{
        .achiral_proof = "fixture demonstrated achirality",
    })).matches);
}

test "coordinate comparisons require an identical run identity" {
    const native = RunIdentity{
        .process_id = "run-4",
        .target = "aarch64-macos",
        .toolchain = "0.17.0-dev.1516+8a4b5424d",
        .optimize_mode = "ReleaseFast",
    };
    try RunIdentity.requireSameBuild(native, native);
    var other_target = native;
    other_target.target = "x86_64-macos";
    try std.testing.expectError(error.MismatchedRunIdentity, RunIdentity.requireSameBuild(native, other_target));
}

test "invariant quality floor rejects a worse native result" {
    const oracle = QualityMetrics{
        .clash_score = 1,
        .bond_length_rms = 0.1,
        .bond_angle_deviation = 0.2,
        .bond_crossings = 2,
        .atoms_inside_rings = 0,
    };
    try std.testing.expect(noWorseThan(oracle, oracle, .{}));
    var worse = oracle;
    worse.bond_crossings += 1;
    try std.testing.expect(!noWorseThan(oracle, worse, .{}));
}
