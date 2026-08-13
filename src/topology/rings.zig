const std = @import("std");
const core = @import("core");
const model = @import("model");

pub const macrocycle_size: u32 = 9;

pub const Flags = packed struct(u8) {
    benzene: bool = false,
    aromatic: bool = false,
    macrocycle: bool = false,
    _padding: u5 = 0,
};

pub const Fusion = struct {
    other: core.ids.RingId,
    atom_start: u32,
    atom_count: u32,
    bond: core.ids.BondId = .invalid,
};

/// Ring chemistry and fusion structure consumed by fragment layout. Ring
/// membership itself is produced during graph preparation; this layer adds the
/// pinned sketcherMinimizerRing/initializeFusedRingInformation semantics.
pub const Analysis = struct {
    allocator: std.mem.Allocator,
    flags: []Flags,
    fusion_offsets: []u32,
    fusions: []Fusion,
    fusion_atoms: []core.ids.AtomId,

    pub fn init(
        allocator: std.mem.Allocator,
        membership: anytype,
        atoms: []const model.Atom,
        bonds: []const model.Bond,
    ) core.errors.Error!Analysis {
        const ring_count = membership.rings.len;
        const flags = allocator.alloc(Flags, ring_count) catch return error.OutOfMemory;
        errdefer allocator.free(flags);
        for (flags, 0..) |*flag, index| {
            const ring_id = core.ids.RingId.fromIndex(@intCast(index));
            flag.* = .{
                .benzene = isBenzene(membership, ring_id, atoms, bonds),
                .aromatic = isAromatic(membership, ring_id, atoms, bonds),
                .macrocycle = membership.atoms(ring_id).len >= macrocycle_size,
            };
        }

        var fusion_list: std.ArrayList(Fusion) = .empty;
        defer fusion_list.deinit(allocator);
        var atom_list: std.ArrayList(core.ids.AtomId) = .empty;
        defer atom_list.deinit(allocator);
        const fusion_offsets = allocator.alloc(u32, ring_count + 1) catch return error.OutOfMemory;
        errdefer allocator.free(fusion_offsets);
        for (0..ring_count) |left_index| {
            fusion_offsets[left_index] = @intCast(fusion_list.items.len);
            const left = core.ids.RingId.fromIndex(@intCast(left_index));
            for (0..ring_count) |right_index| {
                if (left_index == right_index) continue;
                const right = core.ids.RingId.fromIndex(@intCast(right_index));
                const atom_start = atom_list.items.len;
                for (membership.atoms(left)) |atom| {
                    if (std.mem.indexOfScalar(core.ids.AtomId, membership.atoms(right), atom) != null) {
                        atom_list.append(allocator, atom) catch return error.OutOfMemory;
                    }
                }
                const shared_count = atom_list.items.len - atom_start;
                const fusion_bond = if (shared_count == 0)
                    connectingMultipleBond(membership, left, right, bonds)
                else
                    core.ids.BondId.invalid;
                if (shared_count == 0 and fusion_bond == .invalid) continue;
                fusion_list.append(allocator, .{
                    .other = right,
                    .atom_start = @intCast(atom_start),
                    .atom_count = @intCast(shared_count),
                    .bond = fusion_bond,
                }) catch return error.OutOfMemory;
            }
        }
        fusion_offsets[ring_count] = @intCast(fusion_list.items.len);
        const fusions = fusion_list.toOwnedSlice(allocator) catch return error.OutOfMemory;
        errdefer allocator.free(fusions);
        const fusion_atoms = atom_list.toOwnedSlice(allocator) catch return error.OutOfMemory;
        return .{
            .allocator = allocator,
            .flags = flags,
            .fusion_offsets = fusion_offsets,
            .fusions = fusions,
            .fusion_atoms = fusion_atoms,
        };
    }

    pub fn deinit(self: *Analysis) void {
        self.allocator.free(self.fusion_atoms);
        self.allocator.free(self.fusions);
        self.allocator.free(self.fusion_offsets);
        self.allocator.free(self.flags);
        self.* = undefined;
    }

    pub fn fusedWith(self: Analysis, ring: core.ids.RingId) []const Fusion {
        return self.fusions[self.fusion_offsets[ring.index()]..self.fusion_offsets[ring.index() + 1]];
    }

    pub fn fusionAtoms(self: Analysis, fusion: Fusion) []const core.ids.AtomId {
        return self.fusion_atoms[fusion.atom_start..][0..fusion.atom_count];
    }
};

fn isBenzene(membership: anytype, ring: core.ids.RingId, atoms: []const model.Atom, bonds: []const model.Bond) bool {
    const members = membership.atoms(ring);
    if (members.len != 6) return false;
    for (members) |atom| {
        if (atoms[atom.index()].atomic_number != .carbon or !hasDoubleBond(atom, bonds)) return false;
    }
    return true;
}

fn isAromatic(membership: anytype, ring: core.ids.RingId, atoms: []const model.Atom, bonds: []const model.Bond) bool {
    const ring_bonds = membership.ringBonds(ring);
    var double_bonds: u32 = 0;
    for (ring_bonds) |bond| double_bonds += @intFromBool(bonds[bond.index()].effective_order == .double);
    if (ring_bonds.len == 6 and double_bonds == 3) return true;
    if (ring_bonds.len != 5 or double_bonds != 2) return false;
    var nso_count: u32 = 0;
    for (membership.atoms(ring)) |atom| {
        if (hasDoubleBond(atom, bonds)) continue;
        nso_count += @intFromBool(switch (atoms[atom.index()].atomic_number) {
            .nitrogen, .oxygen, .sulfur => true,
            else => false,
        });
    }
    return nso_count == 1;
}

fn hasDoubleBond(atom: core.ids.AtomId, bonds: []const model.Bond) bool {
    for (bonds) |bond| {
        if ((bond.start == atom or bond.end == atom) and bond.effective_order == .double) return true;
    }
    return false;
}

fn connectingMultipleBond(membership: anytype, left: core.ids.RingId, right: core.ids.RingId, bonds: []const model.Bond) core.ids.BondId {
    for (bonds) |bond| {
        if (bond.effective_order == .single or bond.skip) continue;
        const forward = std.mem.indexOfScalar(core.ids.AtomId, membership.atoms(left), bond.start) != null and
            std.mem.indexOfScalar(core.ids.AtomId, membership.atoms(right), bond.end) != null;
        const reverse = std.mem.indexOfScalar(core.ids.AtomId, membership.atoms(left), bond.end) != null and
            std.mem.indexOfScalar(core.ids.AtomId, membership.atoms(right), bond.start) != null;
        if ((forward or reverse) and !shareRing(membership, bond.start, bond.end)) return bond.id;
    }
    return .invalid;
}

fn shareRing(membership: anytype, left: core.ids.AtomId, right: core.ids.AtomId) bool {
    for (membership.atomRings(left)) |ring| {
        if (std.mem.indexOfScalar(core.ids.RingId, membership.atomRings(right), ring) != null) return true;
    }
    return false;
}
