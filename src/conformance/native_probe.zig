//! The native side of the probe surface (cgz-7v2.4.2.1).
//!
//! The epic requires two surfaces around each implementation: a narrow stable
//! one, and a conformance-only one that reports internal state and is never
//! installed. `coordgen_probe.h` is the oracle's; this is native's, and it
//! lives in the conformance layer for the same reason — no production layer may
//! import this module, so nothing installed can reach it.
//!
//! It owns flattened copies of every intermediate seam in the oracle probe:
//! Morgan ranks, rings, fragments, DOFs and penalties, accepted template
//! mappings, components and their placement transforms. It also retains the
//! pose immediately before global orientation for stage localization.

const std = @import("std");
const core = @import("core");
const api = @import("api");
const generator = @import("generator");
const probe_types = @import("probe_types.zig");

const invalid_index = std.math.maxInt(u32);

pub const Probe = struct {
    input_to_internal: []u32,
    internal_to_input: []u32,
    morgan_ranks: []u32,
    ring_atoms: []u32,
    fragment_atoms: []u32,
    fragment_rings: []u32,
    component_atoms: []u32,
    dof_affected_atoms: []u32,
    template_mapping: []probe_types.TemplateMappingProbe,
    rings: []probe_types.RingProbe,
    fragments: []probe_types.FragmentProbe,
    dofs: []probe_types.DofProbe,
    components: []probe_types.ComponentProbe,
    clean_pose: bool,

    fn deinit(self: *Probe, allocator: std.mem.Allocator) void {
        allocator.free(self.components);
        allocator.free(self.dofs);
        allocator.free(self.fragments);
        allocator.free(self.rings);
        allocator.free(self.template_mapping);
        allocator.free(self.dof_affected_atoms);
        allocator.free(self.component_atoms);
        allocator.free(self.fragment_rings);
        allocator.free(self.fragment_atoms);
        allocator.free(self.ring_atoms);
        allocator.free(self.morgan_ranks);
        allocator.free(self.internal_to_input);
        allocator.free(self.input_to_internal);
        self.* = undefined;
    }
};

pub const Stages = struct {
    allocator: std.mem.Allocator,
    /// Final coordinates, caller input order. Identical to what `api.generate`
    /// returns for the same input; the differential runner asserts that rather
    /// than assuming it, which makes any divergence between the two a reported
    /// nondeterminism finding instead of a silent one.
    coordinates: []core.math.Vec2,
    /// The same atoms before `orientComponents` runs.
    pre_orientation: []core.math.Vec2,
    clean_pose: bool,
    /// Counts every component global orientation handled or intentionally
    /// declined. Ring-bearing components belong in `oriented_components`;
    /// only fixed/constrained components belong in the other count.
    orientation: generator.components.OrientationOutcome,
    probe: Probe,

    pub fn deinit(self: *Stages) void {
        self.probe.deinit(self.allocator);
        self.allocator.free(self.pre_orientation);
        self.allocator.free(self.coordinates);
        self.* = undefined;
    }
};

/// Run one generation and report both stages. Validation and the pipeline are
/// the public entry point's, reached through the same `generateInto`, so this
/// cannot drift into measuring a different code path than production runs.
pub fn generateStages(allocator: std.mem.Allocator, input: api.Input) core.errors.Error!Stages {
    try input.validate();

    const coordinates = allocator.alloc(core.math.Vec2, input.atoms.len) catch
        return error.OutOfMemory;
    errdefer allocator.free(coordinates);
    const pre_orientation = allocator.alloc(core.math.Vec2, input.atoms.len) catch
        return error.OutOfMemory;
    errdefer allocator.free(pre_orientation);

    const input_to_internal = allocator.alloc(u32, input.atoms.len) catch
        return error.OutOfMemory;
    defer allocator.free(input_to_internal);
    const internal_to_input = allocator.alloc(core.ids.InputIndex, input.atoms.len) catch
        return error.OutOfMemory;
    defer allocator.free(internal_to_input);

    var orientation: generator.components.OrientationOutcome = .{};
    var capture = Capture{
        .allocator = allocator,
        .pre_orientation = pre_orientation,
    };
    defer capture.deinit();
    const clean_pose = try generator.generateInto(allocator, input, .{
        .coordinates = coordinates,
        .input_to_internal = input_to_internal,
        .internal_to_input = internal_to_input,
        .pre_orientation = pre_orientation,
        .orientation_outcome = &orientation,
        .probe_sink = .{
            .context = &capture,
            .template_observer = .{ .context = &capture, .recordFn = recordTemplate },
            .captureFn = captureSnapshot,
        },
    });
    const probe = capture.probe orelse return error.InvalidMapping;
    capture.probe = null;

    return .{
        .allocator = allocator,
        .coordinates = coordinates,
        .pre_orientation = pre_orientation,
        .clean_pose = clean_pose,
        .orientation = orientation,
        .probe = probe,
    };
}

const TemplateCapture = struct {
    fragment: core.ids.FragmentId,
    template_index: u32,
    mapping_start: u32,
    mapping_count: u32,
};

const Capture = struct {
    allocator: std.mem.Allocator,
    pre_orientation: []const core.math.Vec2,
    template_mapping: std.ArrayList(probe_types.TemplateMappingProbe) = .empty,
    templates: std.ArrayList(TemplateCapture) = .empty,
    probe: ?Probe = null,

    fn deinit(self: *Capture) void {
        if (self.probe) |*probe| probe.deinit(self.allocator);
        self.templates.deinit(self.allocator);
        self.template_mapping.deinit(self.allocator);
    }
};

fn recordTemplate(
    context_ptr: *anyopaque,
    fragment: core.ids.FragmentId,
    template_index: u32,
    atoms: []const core.ids.AtomId,
    mapping: []const u8,
) core.errors.Error!void {
    const capture: *Capture = @ptrCast(@alignCast(context_ptr));
    if (atoms.len != mapping.len or capture.template_mapping.items.len > std.math.maxInt(u32)) {
        return error.InvalidMapping;
    }
    const start: u32 = @intCast(capture.template_mapping.items.len);
    for (atoms, mapping) |atom, template_atom| {
        capture.template_mapping.append(capture.allocator, .{
            // Converted to input order once the canonical map is available.
            .input_atom = atom.index(),
            .template_atom = template_atom,
        }) catch return error.OutOfMemory;
    }
    capture.templates.append(capture.allocator, .{
        .fragment = fragment,
        .template_index = template_index,
        .mapping_start = start,
        .mapping_count = @intCast(mapping.len),
    }) catch return error.OutOfMemory;
}

fn captureSnapshot(context_ptr: *anyopaque, snapshot: generator.ProbeSnapshot) core.errors.Error!void {
    const capture: *Capture = @ptrCast(@alignCast(context_ptr));
    const allocator = capture.allocator;
    const prepared = snapshot.prepared;
    const order = prepared.working.order;

    for (capture.template_mapping.items) |*mapping| {
        if (mapping.input_atom >= order.internal_to_input.len) return error.InvalidMapping;
        mapping.input_atom = order.internal_to_input[mapping.input_atom];
    }

    var ring_atoms: std.ArrayList(u32) = .empty;
    defer ring_atoms.deinit(allocator);
    var rings: std.ArrayList(probe_types.RingProbe) = .empty;
    defer rings.deinit(allocator);
    for (prepared.rings.rings) |ring| {
        const start: u32 = @intCast(ring_atoms.items.len);
        for (prepared.rings.atoms(ring.id)) |atom| {
            ring_atoms.append(allocator, try inputAtom(order.internal_to_input, atom)) catch return error.OutOfMemory;
        }
        rings.append(allocator, .{ .atom_start = start, .atom_count = ring.atom_count }) catch return error.OutOfMemory;
    }

    var fragment_atoms: std.ArrayList(u32) = .empty;
    defer fragment_atoms.deinit(allocator);
    var fragment_rings: std.ArrayList(u32) = .empty;
    defer fragment_rings.deinit(allocator);
    var fragments: std.ArrayList(probe_types.FragmentProbe) = .empty;
    defer fragments.deinit(allocator);
    const fragment_order = allocator.alloc(core.ids.FragmentId, snapshot.fragmentation.fragments.len) catch return error.OutOfMemory;
    defer allocator.free(fragment_order);
    const fragment_probe_id = allocator.alloc(u32, snapshot.fragmentation.fragments.len) catch return error.OutOfMemory;
    defer allocator.free(fragment_probe_id);
    @memset(fragment_probe_id, invalid_index);
    var fragment_order_len: usize = 0;
    for (0..prepared.graph.component_count) |component_index| {
        const component = core.ids.MoleculeId.fromIndex(@intCast(component_index));
        const root = snapshot.fragmentation.main_fragments[component.index()];
        fragment_order[fragment_order_len] = root;
        fragment_probe_id[root.index()] = @intCast(fragment_order_len);
        fragment_order_len += 1;
        var queue_index = fragment_order_len - 1;
        while (queue_index < fragment_order_len) : (queue_index += 1) {
            const parent = fragment_order[queue_index];
            for (prepared.graph.structuralBonds()) |bond_id| {
                for (snapshot.fragmentation.fragments) |candidate| {
                    if (candidate.parent != parent or fragment_probe_id[candidate.id.index()] != invalid_index) continue;
                    if (candidate.bond_to_parent != bond_id) continue;
                    fragment_order[fragment_order_len] = candidate.id;
                    fragment_probe_id[candidate.id.index()] = @intCast(fragment_order_len);
                    fragment_order_len += 1;
                }
            }
        }
    }
    if (fragment_order_len != snapshot.fragmentation.fragments.len) return error.InvalidMapping;
    const emitted_atoms = allocator.alloc(bool, prepared.working.atoms.len) catch return error.OutOfMemory;
    defer allocator.free(emitted_atoms);
    for (fragment_order) |fragment_id| {
        const fragment = snapshot.fragmentation.fragments[fragment_id.index()];
        const atom_start: u32 = @intCast(fragment_atoms.items.len);
        @memset(emitted_atoms, false);
        for (prepared.graph.structuralBonds()) |bond_id| {
            const bond = prepared.working.bonds[bond_id.index()];
            for ([_]core.ids.AtomId{ bond.start, bond.end }) |atom| {
                if (snapshot.fragmentation.atom_fragment[atom.index()] != fragment.id or emitted_atoms[atom.index()]) continue;
                fragment_atoms.append(allocator, try inputAtom(order.internal_to_input, atom)) catch return error.OutOfMemory;
                emitted_atoms[atom.index()] = true;
            }
        }
        for (prepared.graph.componentMembers(fragment.component)) |atom| {
            if (snapshot.fragmentation.atom_fragment[atom.index()] == fragment.id and !emitted_atoms[atom.index()]) {
                fragment_atoms.append(allocator, try inputAtom(order.internal_to_input, atom)) catch return error.OutOfMemory;
            }
        }
        const ring_start: u32 = @intCast(fragment_rings.items.len);
        for (prepared.rings.rings) |ring| {
            const members = prepared.rings.atoms(ring.id);
            if (members.len != 0 and snapshot.fragmentation.atom_fragment[members[0].index()] == fragment.id) {
                fragment_rings.append(allocator, ring.id.index()) catch return error.OutOfMemory;
            }
        }
        var dof_start: u32 = 0;
        var dof_count: u32 = 0;
        for (fragment_order[0..fragment_probe_id[fragment.id.index()]]) |prior_fragment| {
            for (snapshot.dofs.items) |dof| {
                if (dof.fragment == prior_fragment) dof_start += 1;
            }
        }
        for (snapshot.dofs.items) |dof| {
            if (dof.fragment == fragment.id) dof_count += 1;
        }
        var flags: u32 = @intFromBool(fragment.flags.fixed) |
            (@as(u32, @intFromBool(fragment.flags.constrained)) << 2) |
            (@as(u32, @intFromBool(fragment.flags.constrained_flip)) << 3) |
            (@as(u32, @intFromBool(fragment.flags.chain)) << 4);
        var template_match: u32 = invalid_index;
        var mapping_start: u32 = 0;
        var mapping_count: u32 = 0;
        for (capture.templates.items) |template| if (template.fragment == fragment.id) {
            flags |= 1 << 1;
            template_match = template.template_index;
            mapping_start = template.mapping_start;
            mapping_count = template.mapping_count;
            break;
        };
        fragments.append(allocator, .{
            .parent = if (fragment.parent.isValid()) fragment_probe_id[fragment.parent.index()] else invalid_index,
            .component = fragment.component.index(),
            .atom_start = atom_start,
            .atom_count = @intCast(fragment_atoms.items.len - atom_start),
            .ring_start = ring_start,
            .ring_count = @intCast(fragment_rings.items.len - ring_start),
            .dof_start = dof_start,
            .dof_count = dof_count,
            .flags = flags,
            .template_match = template_match,
            .template_mapping_start = mapping_start,
            .template_mapping_count = mapping_count,
        }) catch return error.OutOfMemory;
    }

    var affected: std.ArrayList(u32) = .empty;
    defer affected.deinit(allocator);
    var dofs: std.ArrayList(probe_types.DofProbe) = .empty;
    defer dofs.deinit(allocator);
    for (fragment_order) |fragment_id| {
        for (snapshot.dofs.items) |dof| {
            if (dof.fragment != fragment_id) continue;
            const affected_start: u32 = @intCast(affected.items.len);
            const range = snapshot.dofs.affected_atoms[dof.affected_atoms.start..][0..dof.affected_atoms.len];
            for (range) |atom| affected.append(allocator, try inputAtom(order.internal_to_input, atom)) catch return error.OutOfMemory;
            const fragment = snapshot.fragmentation.fragments[dof.fragment.index()];
            const chain_parent = fragment.parent.isValid() and fragment.flags.chain and
                snapshot.fragmentation.fragments[fragment.parent.index()].flags.chain;
            const penalty = @import("optimize").discrete.penalty(dof, range.len, .{
                .constrained_flip = fragment.flags.constrained_flip,
                .chain_with_chain_parent = chain_parent,
            }) catch return error.InvalidMapping;
            var probe = probe_types.dofProbe(dof, penalty);
            probe.id = @intCast(dofs.items.len);
            probe.fragment = fragment_probe_id[dof.fragment.index()];
            probe.affected_start = affected_start;
            probe.atom_a = try mapOptionalAtom(order.internal_to_input, probe.atom_a);
            probe.atom_b = try mapOptionalAtom(order.internal_to_input, probe.atom_b);
            dofs.append(allocator, probe) catch return error.OutOfMemory;
        }
    }

    var component_atoms: std.ArrayList(u32) = .empty;
    defer component_atoms.deinit(allocator);
    var components: std.ArrayList(probe_types.ComponentProbe) = .empty;
    defer components.deinit(allocator);
    for (0..prepared.graph.component_count) |index| {
        const members = prepared.graph.componentMembers(core.ids.MoleculeId.fromIndex(@intCast(index)));
        const start: u32 = @intCast(component_atoms.items.len);
        for (members) |atom| component_atoms.append(allocator, try inputAtom(order.internal_to_input, atom)) catch return error.OutOfMemory;
        components.append(allocator, .{
            .atom_start = start,
            .atom_count = @intCast(members.len),
            .transform_status = .observed,
            .reserved = 0,
            .transform = try fitTransform(capture.pre_orientation, prepared.working.atoms, order.internal_to_input, members),
        }) catch return error.OutOfMemory;
    }

    const owned_input_to_internal = allocator.dupe(u32, order.input_to_internal) catch return error.OutOfMemory;
    errdefer allocator.free(owned_input_to_internal);
    const owned_internal_to_input = try copyInputIndices(allocator, order.internal_to_input);
    errdefer allocator.free(owned_internal_to_input);
    const owned_morgan_ranks = allocator.dupe(u32, prepared.component_morgan_scores) catch return error.OutOfMemory;
    errdefer allocator.free(owned_morgan_ranks);
    const owned_ring_atoms = ring_atoms.toOwnedSlice(allocator) catch return error.OutOfMemory;
    errdefer allocator.free(owned_ring_atoms);
    const owned_fragment_atoms = fragment_atoms.toOwnedSlice(allocator) catch return error.OutOfMemory;
    errdefer allocator.free(owned_fragment_atoms);
    const owned_fragment_rings = fragment_rings.toOwnedSlice(allocator) catch return error.OutOfMemory;
    errdefer allocator.free(owned_fragment_rings);
    const owned_component_atoms = component_atoms.toOwnedSlice(allocator) catch return error.OutOfMemory;
    errdefer allocator.free(owned_component_atoms);
    const owned_affected = affected.toOwnedSlice(allocator) catch return error.OutOfMemory;
    errdefer allocator.free(owned_affected);
    const owned_template_mapping = capture.template_mapping.toOwnedSlice(allocator) catch return error.OutOfMemory;
    errdefer allocator.free(owned_template_mapping);
    const owned_rings = rings.toOwnedSlice(allocator) catch return error.OutOfMemory;
    errdefer allocator.free(owned_rings);
    const owned_fragments = fragments.toOwnedSlice(allocator) catch return error.OutOfMemory;
    errdefer allocator.free(owned_fragments);
    const owned_dofs = dofs.toOwnedSlice(allocator) catch return error.OutOfMemory;
    errdefer allocator.free(owned_dofs);
    const owned_components = components.toOwnedSlice(allocator) catch return error.OutOfMemory;
    errdefer allocator.free(owned_components);
    capture.probe = .{
        .input_to_internal = owned_input_to_internal,
        .internal_to_input = owned_internal_to_input,
        .morgan_ranks = owned_morgan_ranks,
        .ring_atoms = owned_ring_atoms,
        .fragment_atoms = owned_fragment_atoms,
        .fragment_rings = owned_fragment_rings,
        .component_atoms = owned_component_atoms,
        .dof_affected_atoms = owned_affected,
        .template_mapping = owned_template_mapping,
        .rings = owned_rings,
        .fragments = owned_fragments,
        .dofs = owned_dofs,
        .components = owned_components,
        .clean_pose = snapshot.clean_pose,
    };
}

fn inputAtom(map: []const core.ids.InputIndex, atom: core.ids.AtomId) core.errors.Error!u32 {
    if (!atom.isValid() or atom.index() >= map.len) return error.InvalidMapping;
    return map[atom.index()];
}

fn mapOptionalAtom(map: []const core.ids.InputIndex, atom: u32) core.errors.Error!u32 {
    if (atom == invalid_index) return invalid_index;
    if (atom >= map.len) return error.InvalidMapping;
    return map[atom];
}

fn copyInputIndices(allocator: std.mem.Allocator, values: []const core.ids.InputIndex) core.errors.Error![]u32 {
    const output = allocator.alloc(u32, values.len) catch return error.OutOfMemory;
    for (values, output) |value, *slot| slot.* = value;
    return output;
}

fn fitTransform(
    before: []const core.math.Vec2,
    after: []const @import("model").Atom,
    map: []const core.ids.InputIndex,
    members: []const core.ids.AtomId,
) core.errors.Error![6]f32 {
    if (members.len == 0) return error.InvalidMapping;
    var source_center: core.math.Vec2 = .{};
    var target_center: core.math.Vec2 = .{};
    for (members) |atom| {
        const input_index = try inputAtom(map, atom);
        source_center.x += before[input_index].x;
        source_center.y += before[input_index].y;
        target_center.x += after[atom.index()].coordinates.x;
        target_center.y += after[atom.index()].coordinates.y;
    }
    const scale = 1 / @as(f32, @floatFromInt(members.len));
    source_center.x *= scale;
    source_center.y *= scale;
    target_center.x *= scale;
    target_center.y *= scale;
    var dot: f32 = 0;
    var cross: f32 = 0;
    var reflected_dot: f32 = 0;
    var reflected_cross: f32 = 0;
    for (members) |atom| {
        const input_index = try inputAtom(map, atom);
        const sx = before[input_index].x - source_center.x;
        const sy = before[input_index].y - source_center.y;
        const tx = after[atom.index()].coordinates.x - target_center.x;
        const ty = after[atom.index()].coordinates.y - target_center.y;
        dot += sx * tx + sy * ty;
        cross += sx * ty - sy * tx;
        reflected_dot += sx * tx - sy * ty;
        reflected_cross += sy * tx + sx * ty;
    }
    const norm = @sqrt(dot * dot + cross * cross);
    const cosine: f32 = if (norm == 0) 1 else dot / norm;
    const sine: f32 = if (norm == 0) 0 else cross / norm;
    const direct = [6]f32{
        cosine,
        -sine,
        sine,
        cosine,
        target_center.x - cosine * source_center.x + sine * source_center.y,
        target_center.y - sine * source_center.x - cosine * source_center.y,
    };
    const reflected_norm = @sqrt(reflected_dot * reflected_dot + reflected_cross * reflected_cross);
    const reflected_cosine: f32 = if (reflected_norm == 0) 1 else reflected_dot / reflected_norm;
    const reflected_sine: f32 = if (reflected_norm == 0) 0 else reflected_cross / reflected_norm;
    const reflected = [6]f32{
        reflected_cosine,
        reflected_sine,
        reflected_sine,
        -reflected_cosine,
        target_center.x - reflected_cosine * source_center.x - reflected_sine * source_center.y,
        target_center.y - reflected_sine * source_center.x + reflected_cosine * source_center.y,
    };
    return if (transformResidual(before, after, map, members, reflected) <
        transformResidual(before, after, map, members, direct)) reflected else direct;
}

fn transformResidual(
    before: []const core.math.Vec2,
    after: []const @import("model").Atom,
    map: []const core.ids.InputIndex,
    members: []const core.ids.AtomId,
    transform: [6]f32,
) f32 {
    var residual: f32 = 0;
    for (members) |atom| {
        const input_index = map[atom.index()];
        const source = before[input_index];
        const target = after[atom.index()].coordinates;
        const dx = transform[0] * source.x + transform[1] * source.y + transform[4] - target.x;
        const dy = transform[2] * source.x + transform[3] * source.y + transform[5] - target.y;
        residual += dx * dx + dy * dy;
    }
    return residual;
}

test "the probe reports a pre-orientation pose that global orientation then moves" {
    const atoms = [_]api.AtomInput{
        .{ .atomic_number = .carbon },
        .{ .atomic_number = .carbon },
        .{ .atomic_number = .carbon },
        .{ .atomic_number = .carbon },
        .{ .atomic_number = .carbon },
        .{ .atomic_number = .carbon },
        .{ .atomic_number = .oxygen },
    };
    const bonds = [_]api.BondInput{
        .{ .start = 0, .end = 1, .order = .single },
        .{ .start = 1, .end = 2, .order = .single },
        .{ .start = 2, .end = 3, .order = .single },
        .{ .start = 3, .end = 4, .order = .single },
        .{ .start = 4, .end = 5, .order = .single },
        .{ .start = 5, .end = 0, .order = .single },
        .{ .start = 2, .end = 6, .order = .single },
    };
    const input = api.Input{ .atoms = &atoms, .bonds = &bonds };

    var stages = try generateStages(std.testing.allocator, input);
    defer stages.deinit();

    // This component contains a six-membered ring. The count is the runtime
    // declaration that prevents the old `has_ring` decline from returning
    // unnoticed: the only legal decline category is fixed/constrained.
    try std.testing.expectEqual(@as(u32, 1), stages.orientation.oriented_components);
    try std.testing.expectEqual(@as(u32, 0), stages.orientation.constrained_components);
    try std.testing.expectEqual(atoms.len, stages.pre_orientation.len);
    for (stages.pre_orientation) |point| try std.testing.expect(point.isFinite());
    for (stages.coordinates) |point| try std.testing.expect(point.isFinite());

    // Global orientation is rigid — rotation, reflection and translation only —
    // so it cannot change any distance within a component. Comparing the two
    // stages pairwise is therefore both the strongest available invariant and
    // the assertion that separates a real captured pose from zeroed storage,
    // which would pass a finiteness check.
    for (0..atoms.len) |i| {
        for (i + 1..atoms.len) |j| {
            try std.testing.expectApproxEqAbs(
                separation(stages.pre_orientation, i, j),
                separation(stages.coordinates, i, j),
                0.01,
            );
        }
    }
    // The captured pose is a real layout, not the origin: bonded atoms sit
    // about one bond length apart in it.
    for (bonds) |bond| {
        const length = separation(stages.pre_orientation, bond.start, bond.end);
        try std.testing.expect(length > core.math.bond_length * 0.5);
        try std.testing.expect(length < core.math.bond_length * 1.5);
    }

    // And the final pose agrees with the public entry point's, which is the
    // invariant that lets the runner trust one generation for both.
    var public = try api.generate(std.testing.allocator, input);
    defer public.deinit();
    try std.testing.expectEqual(public.clean_pose, stages.clean_pose);
    for (public.coordinates, stages.coordinates) |expected, actual| {
        try std.testing.expectEqual(expected.x, actual.x);
        try std.testing.expectEqual(expected.y, actual.y);
    }
}

fn separation(coordinates: []const core.math.Vec2, left: usize, right: usize) f32 {
    const delta_x = coordinates[left].x - coordinates[right].x;
    const delta_y = coordinates[left].y - coordinates[right].y;
    return @sqrt(delta_x * delta_x + delta_y * delta_y);
}
