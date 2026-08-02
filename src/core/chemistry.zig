const std = @import("std");

/// Atomic number 0 is reserved for internal virtual atoms. Public molecule
/// input accepts the IUPAC range 1...118 only.
pub const AtomicNumber = enum(u32) {
    virtual = 0,
    hydrogen = 1,
    helium = 2,
    lithium = 3,
    beryllium = 4,
    boron = 5,
    carbon = 6,
    nitrogen = 7,
    oxygen = 8,
    fluorine = 9,
    neon = 10,
    sodium = 11,
    magnesium = 12,
    aluminium = 13,
    silicon = 14,
    phosphorus = 15,
    sulfur = 16,
    chlorine = 17,
    argon = 18,
    potassium = 19,
    calcium = 20,
    scandium = 21,
    titanium = 22,
    vanadium = 23,
    chromium = 24,
    manganese = 25,
    iron = 26,
    cobalt = 27,
    nickel = 28,
    copper = 29,
    zinc = 30,
    gallium = 31,
    germanium = 32,
    arsenic = 33,
    selenium = 34,
    bromine = 35,
    krypton = 36,
    rubidium = 37,
    strontium = 38,
    yttrium = 39,
    zirconium = 40,
    niobium = 41,
    molybdenum = 42,
    technetium = 43,
    ruthenium = 44,
    rhodium = 45,
    palladium = 46,
    silver = 47,
    cadmium = 48,
    indium = 49,
    tin = 50,
    antimony = 51,
    tellurium = 52,
    iodine = 53,
    xenon = 54,
    caesium = 55,
    barium = 56,
    lanthanum = 57,
    cerium = 58,
    praseodymium = 59,
    neodymium = 60,
    promethium = 61,
    samarium = 62,
    europium = 63,
    gadolinium = 64,
    terbium = 65,
    dysprosium = 66,
    holmium = 67,
    erbium = 68,
    thulium = 69,
    ytterbium = 70,
    lutetium = 71,
    hafnium = 72,
    tantalum = 73,
    tungsten = 74,
    rhenium = 75,
    osmium = 76,
    iridium = 77,
    platinum = 78,
    gold = 79,
    mercury = 80,
    thallium = 81,
    lead = 82,
    bismuth = 83,
    polonium = 84,
    astatine = 85,
    radon = 86,
    francium = 87,
    radium = 88,
    actinium = 89,
    thorium = 90,
    protactinium = 91,
    uranium = 92,
    neptunium = 93,
    plutonium = 94,
    americium = 95,
    curium = 96,
    berkelium = 97,
    californium = 98,
    einsteinium = 99,
    fermium = 100,
    mendelevium = 101,
    nobelium = 102,
    lawrencium = 103,
    rutherfordium = 104,
    dubnium = 105,
    seaborgium = 106,
    bohrium = 107,
    hassium = 108,
    meitnerium = 109,
    darmstadtium = 110,
    roentgenium = 111,
    copernicium = 112,
    nihonium = 113,
    flerovium = 114,
    moscovium = 115,
    livermorium = 116,
    tennessine = 117,
    oganesson = 118,

    pub fn fromPublic(value: u32) ?AtomicNumber {
        if (value < 1 or value > 118) return null;
        return @fromBackingInt(@intCast(value));
    }
};

pub const FormalCharge = i32;

pub const BondOrder = enum(u32) {
    zero = 0,
    single = 1,
    double = 2,
    triple = 3,

    pub fn fromInt(value: u32) ?BondOrder {
        return switch (value) {
            0 => .zero,
            1 => .single,
            2 => .double,
            3 => .triple,
            else => null,
        };
    }
};

/// Both upstream relative input and its derived absolute label are preserved.
pub const AtomStereo = enum(u32) {
    unspecified = 0,
    clockwise = 1,
    counter_clockwise = 2,
    r = 3,
    s = 4,
};

/// `cis`/`trans` describe referenced substituents; `z`/`e` describe the
/// absolute double-bond label after CIP conversion.
pub const BondStereo = enum(u32) {
    unspecified = 0,
    cis = 1,
    trans = 2,
    z = 3,
    e = 4,
};

/// Direction is part of the display value so the three upstream booleans
/// (has display, wedge kind, reversed) cannot form inconsistent states.
pub const BondDisplay = enum(u32) {
    none = 0,
    solid_forward = 1,
    solid_reverse = 2,
    hashed_forward = 3,
    hashed_reverse = 4,
};

test "chemical enums have stable widths and zero-order bonds are conserved" {
    try std.testing.expectEqual(@as(usize, 4), @sizeOf(AtomicNumber));
    try std.testing.expectEqual(@as(usize, 4), @sizeOf(BondOrder));
    try std.testing.expectEqual(@as(u32, 0), @backingInt(BondOrder.zero));
    try std.testing.expectEqual(AtomicNumber.oganesson, AtomicNumber.fromPublic(118).?);
    try std.testing.expect(AtomicNumber.fromPublic(0) == null);
    try std.testing.expect(BondOrder.fromInt(4) == null);
}
