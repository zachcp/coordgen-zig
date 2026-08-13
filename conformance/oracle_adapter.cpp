#include "coordgen_abi.h"
#include "coordgen_probe.h"
#include "coordgen_oracle_hook.hpp"

#include "sketcherMinimizer.h"

#include <cmath>
#include <cstdint>
#include <cstring>
#include <new>
#include <unordered_map>
#include <vector>

namespace {

struct Generation {
    std::vector<coordgen_vec2_t> coordinates;
    std::vector<uint32_t> input_to_internal;
    std::vector<uint32_t> internal_to_input;
    std::vector<uint32_t> effective_bond_orders;
    std::vector<uint32_t> bond_displays;
    std::vector<uint32_t> atom_stereo;
    std::vector<uint32_t> morgan_ranks;
    std::vector<uint32_t> ring_atoms;
    std::vector<uint32_t> fragment_atoms;
    std::vector<uint32_t> fragment_rings;
    std::vector<uint32_t> component_atoms;
    std::vector<coordgen_probe_template_mapping_t> template_mapping;
    std::vector<coordgen_probe_ring_t> rings;
    std::vector<coordgen_probe_fragment_t> fragments;
    std::vector<coordgen_probe_dof_t> dofs;
    std::vector<coordgen_probe_component_t> components;
    uint32_t clean_pose = 0;
};

struct TemplateCapture {
    sketcherMinimizerFragment* fragment = nullptr;
    uint32_t template_index = COORDGEN_INVALID_INDEX;
    std::vector<coordgen_probe_template_mapping_t> mapping;
};

struct ComponentCoordinates {
    sketcherMinimizerMolecule* component = nullptr;
    std::vector<coordgen_vec2_t> coordinates;
};

class OracleCapture final : public coordgen_oracle_hook::CaptureSink {
  public:
    std::vector<TemplateCapture> templates;
    std::vector<ComponentCoordinates> before;
    std::vector<ComponentCoordinates> after;

    void recordTemplateMatch(sketcherMinimizerFragment* fragment,
                             const std::vector<sketcherMinimizerAtom*>& atoms,
                             std::size_t template_index,
                             const std::vector<unsigned int>& mapping) override
    {
        TemplateCapture capture;
        capture.fragment = fragment;
        capture.template_index = static_cast<uint32_t>(template_index);
        capture.mapping.reserve(mapping.size());
        for (std::size_t index = 0; index < mapping.size(); ++index) {
            const sketcherMinimizerAtom* atom = atoms[index];
            const uint32_t input_atom = atom->m_chmN < 0
                                           ? COORDGEN_INVALID_INDEX
                                           : static_cast<uint32_t>(atom->m_chmN);
            capture.mapping.push_back({ input_atom, mapping[index] });
        }
        templates.push_back(std::move(capture));
    }

    void captureComponentTransforms(
        const std::vector<sketcherMinimizerMolecule*>& components,
        bool after_placement) override
    {
        std::vector<ComponentCoordinates>& target = after_placement ? after : before;
        target.clear();
        target.reserve(components.size());
        for (sketcherMinimizerMolecule* component : components) {
            ComponentCoordinates capture;
            capture.component = component;
            capture.coordinates.reserve(component->_atoms.size());
            for (const sketcherMinimizerAtom* atom : component->_atoms) {
                const sketcherMinimizerPointF& point = atom->getCoordinates();
                capture.coordinates.push_back({ point.x(), point.y() });
            }
            target.push_back(std::move(capture));
        }
    }
};

struct CaptureScope {
    explicit CaptureScope(OracleCapture* capture) { coordgen_oracle_hook::setCaptureSink(capture); }
    ~CaptureScope() { coordgen_oracle_hook::setCaptureSink(nullptr); }
};

/* initialize() adopts the molecule. Until that call, upstream's molecule
 * destructor does not own the atoms/bonds it contains, so the adapter supplies
 * the missing exception-safe cleanup at the boundary. */
struct PendingMolecule {
    sketcherMinimizerMolecule* value = new sketcherMinimizerMolecule();

    ~PendingMolecule() {
        if (value == nullptr) return;
        for (sketcherMinimizerAtom* atom : value->_atoms) delete atom;
        for (sketcherMinimizerBond* bond : value->_bonds) delete bond;
        delete value;
    }

    sketcherMinimizerMolecule* operator->() const { return value; }
    sketcherMinimizerMolecule* release() {
        sketcherMinimizerMolecule* released = value;
        value = nullptr;
        return released;
    }
};

bool validFlag(uint32_t value) { return value == 0 || value == 1; }

bool finite(float value) { return std::isfinite(value); }

coordgen_error_t validateInput(const coordgen_input_t* input) {
    if (input == nullptr || input->atoms.len == 0) return COORDGEN_ERROR_EMPTY_GRAPH;
    if (input->atoms.ptr == nullptr) return COORDGEN_ERROR_EMPTY_GRAPH;
    if ((input->bonds.len != 0 && input->bonds.ptr == nullptr) ||
        (input->residues.len != 0 && input->residues.ptr == nullptr) ||
        (input->residue_interactions.len != 0 && input->residue_interactions.ptr == nullptr) ||
        (input->extra_bonds.len != 0 && input->extra_bonds.ptr == nullptr)) {
        return COORDGEN_ERROR_INVALID_MAPPING;
    }
    if (input->residues.len != 0 || input->residue_interactions.len != 0 ||
        input->extra_bonds.len != 0) {
        return COORDGEN_ERROR_UNSUPPORTED;
    }
    if (!finite(input->options.precision) || input->options.precision <= 0.0f ||
        !validFlag(input->options.score_residue_interactions) ||
        !validFlag(input->options.treat_nonterminal_bonds_to_metal_as_zero_order) ||
        !validFlag(input->options.even_angles) ||
        !validFlag(input->options.skip_minimization) ||
        !validFlag(input->options.force_open_macrocycles) ||
        !validFlag(input->options.constrain_all_atoms) ||
        !validFlag(input->options.build_from_fragments) ||
        !validFlag(input->options.debug_coordinates) ||
        !validFlag(input->options.load_templates)) {
        return COORDGEN_ERROR_INVALID_OPTION;
    }
    if (input->options.template_directory.len != 0 || input->options.template_directory.ptr != nullptr) {
        /* A template directory is process-global upstream state.  The adapter
         * deliberately refuses it until its lifetime and concurrency contract
         * can be made explicit. */
        return COORDGEN_ERROR_UNSUPPORTED;
    }
    if (input->options.build_from_fragments != 0) {
        /* Upstream's sketcherMinimizer::buildFromFragments(bool) is not a
         * stored option: it forwards straight to
         * CoordgenMinimizer::buildFromFragments(bool firstTime) const, an
         * imperative pipeline step, not a flag any later stage reads back.
         * Every real upstream call site (sketcherMinimizer.cpp:300, 1160,
         * 2287) invokes it with firstTime=true from inside a sequence
         * (findFragments(); buildFromFragments(true); avoidClashes(); ...)
         * that runGenerateCoordinates() already performs unconditionally.
         * There is no upstream conditional this ABI flag could faithfully
         * gate within the single-shot coordgen_generate() path: accepting
         * true and silently no-op'ing it (as the previous pre-initialize()
         * call did) misrepresents the option as doing something it cannot
         * do without inventing pipeline behavior upstream itself never
         * runs. The default (0/false) is accepted because it matches
         * actual pinned behavior: fragments are always built during
         * generation regardless of this flag. See cgz-7v2.8. */
        return COORDGEN_ERROR_UNSUPPORTED;
    }
    for (uint32_t i = 0; i < input->atoms.len; ++i) {
        const coordgen_atom_input_t& atom = input->atoms.ptr[i];
        if (atom.atomic_number == 0 || atom.atomic_number > 118) return COORDGEN_ERROR_INVALID_ATOMIC_NUMBER;
        if (atom.stereo > COORDGEN_ATOM_STEREO_S) return COORDGEN_ERROR_INVALID_STEREO;
        if ((atom.flags & COORDGEN_ATOM_HAS_TEMPLATE_COORDINATES) != 0 &&
            (!finite(atom.template_coordinates.x) || !finite(atom.template_coordinates.y))) {
            return COORDGEN_ERROR_INVALID_COORDINATE;
        }
        if ((atom.flags & COORDGEN_ATOM_HAS_COORDINATES_3D) != 0 &&
            (!finite(atom.coordinates_3d.x) || !finite(atom.coordinates_3d.y) || !finite(atom.coordinates_3d.z))) {
            return COORDGEN_ERROR_INVALID_COORDINATE;
        }
        if (atom.stereo != COORDGEN_ATOM_STEREO_UNSPECIFIED &&
            ((atom.stereo_looking_from != COORDGEN_INVALID_INDEX && atom.stereo_looking_from >= input->atoms.len) ||
             (atom.stereo_atom_a != COORDGEN_INVALID_INDEX && atom.stereo_atom_a >= input->atoms.len) ||
             (atom.stereo_atom_b != COORDGEN_INVALID_INDEX && atom.stereo_atom_b >= input->atoms.len))) {
            return COORDGEN_ERROR_INVALID_ATOM_INDEX;
        }
    }
    for (uint32_t i = 0; i < input->bonds.len; ++i) {
        const coordgen_bond_input_t& bond = input->bonds.ptr[i];
        if (bond.start >= input->atoms.len || bond.end >= input->atoms.len) return COORDGEN_ERROR_INVALID_ATOM_INDEX;
        if (bond.order > COORDGEN_BOND_TRIPLE) return COORDGEN_ERROR_INVALID_BOND_ORDER;
        if (bond.stereo > COORDGEN_BOND_STEREO_E || bond.display > COORDGEN_BOND_DISPLAY_HASHED_REVERSE) {
            return COORDGEN_ERROR_INVALID_STEREO;
        }
        if (!finite(bond.crossing_penalty_multiplier)) return COORDGEN_ERROR_INVALID_COORDINATE;
        if (bond.stereo != COORDGEN_BOND_STEREO_UNSPECIFIED &&
            ((bond.stereo_atom_a != COORDGEN_INVALID_INDEX && bond.stereo_atom_a >= input->atoms.len) ||
             (bond.stereo_atom_b != COORDGEN_INVALID_INDEX && bond.stereo_atom_b >= input->atoms.len))) {
            return COORDGEN_ERROR_INVALID_ATOM_INDEX;
        }
    }
    return COORDGEN_OK;
}

uint32_t atomIndex(const sketcherMinimizerAtom* atom) {
    return atom != nullptr && atom->m_chmN >= 0 ? static_cast<uint32_t>(atom->m_chmN) : COORDGEN_INVALID_INDEX;
}

uint32_t bondDisplay(const sketcherMinimizerBond* bond) {
    if (!bond->hasStereochemistryDisplay) return COORDGEN_BOND_DISPLAY_NONE;
    if (bond->isWedge) return bond->isReversed ? COORDGEN_BOND_DISPLAY_SOLID_REVERSE : COORDGEN_BOND_DISPLAY_SOLID_FORWARD;
    return bond->isReversed ? COORDGEN_BOND_DISPLAY_HASHED_REVERSE : COORDGEN_BOND_DISPLAY_HASHED_FORWARD;
}

uint32_t atomStereo(const sketcherMinimizerAtom* atom) {
    if (!atom->m_isStereogenic) return COORDGEN_ATOM_STEREO_UNSPECIFIED;
    return atom->isR ? COORDGEN_ATOM_STEREO_R : COORDGEN_ATOM_STEREO_S;
}

uint32_t dofKind(const CoordgenFragmentDOF* dof) {
    if (dynamic_cast<const CoordgenFlipFragmentDOF*>(dof) != nullptr) return COORDGEN_PROBE_DOF_FLIP_FRAGMENT;
    if (dynamic_cast<const CoordgenChangeParentBondLengthFragmentDOF*>(dof) != nullptr) return COORDGEN_PROBE_DOF_CHANGE_PARENT_BOND_LENGTH;
    if (dynamic_cast<const CoordgenRotateFragmentDOF*>(dof) != nullptr) return COORDGEN_PROBE_DOF_ROTATE_FRAGMENT;
    if (dynamic_cast<const CoordgenScaleAtomsDOF*>(dof) != nullptr) return COORDGEN_PROBE_DOF_SCALE_ATOMS;
    if (dynamic_cast<const CoordgenScaleFragmentDOF*>(dof) != nullptr) return COORDGEN_PROBE_DOF_SCALE_FRAGMENT;
    if (dynamic_cast<const CoordgenInvertBondDOF*>(dof) != nullptr) return COORDGEN_PROBE_DOF_INVERT_BOND;
    return COORDGEN_PROBE_DOF_FLIP_RING;
}

void setObservedTransform(const ComponentCoordinates* before,
                          const ComponentCoordinates* after,
                          coordgen_probe_component_t& result)
{
    if (before == nullptr || after == nullptr || before->coordinates.empty() ||
        before->coordinates.size() != after->coordinates.size()) {
        return;
    }
    result.transform_status = COORDGEN_PROBE_TRANSFORM_OBSERVED;
    float before_x = 0.0f;
    float before_y = 0.0f;
    float after_x = 0.0f;
    float after_y = 0.0f;
    for (std::size_t i = 0; i < before->coordinates.size(); ++i) {
        before_x += before->coordinates[i].x;
        before_y += before->coordinates[i].y;
        after_x += after->coordinates[i].x;
        after_y += after->coordinates[i].y;
    }
    const float count = static_cast<float>(before->coordinates.size());
    before_x /= count;
    before_y /= count;
    after_x /= count;
    after_y /= count;
    float xx = 0.0f;
    float xy = 0.0f;
    float yy = 0.0f;
    float ax = 0.0f;
    float ay = 0.0f;
    float bx = 0.0f;
    float by = 0.0f;
    for (std::size_t i = 0; i < before->coordinates.size(); ++i) {
        const float x = before->coordinates[i].x - before_x;
        const float y = before->coordinates[i].y - before_y;
        const float u = after->coordinates[i].x - after_x;
        const float v = after->coordinates[i].y - after_y;
        xx += x * x;
        xy += x * y;
        yy += y * y;
        ax += u * x;
        ay += u * y;
        bx += v * x;
        by += v * y;
    }
    float m00 = 1.0f;
    float m01 = 0.0f;
    float m10 = 0.0f;
    float m11 = 1.0f;
    const float determinant = xx * yy - xy * xy;
    if (std::fabs(determinant) > 0.000001f) {
        m00 = (ax * yy - ay * xy) / determinant;
        m01 = (ay * xx - ax * xy) / determinant;
        m10 = (bx * yy - by * xy) / determinant;
        m11 = (by * xx - bx * xy) / determinant;
    }
    result.transform[0] = m00;
    result.transform[1] = m01;
    result.transform[2] = m10;
    result.transform[3] = m11;
    result.transform[4] = after_x - (m00 * before_x + m01 * before_y);
    result.transform[5] = after_y - (m10 * before_x + m11 * before_y);
}

void fillStable(Generation& value, coordgen_result_t* result) {
    result->coordinates = { value.coordinates.data(), static_cast<uint32_t>(value.coordinates.size()), 0 };
    result->input_to_internal = { value.input_to_internal.data(), static_cast<uint32_t>(value.input_to_internal.size()), 0 };
    result->internal_to_input = { value.internal_to_input.data(), static_cast<uint32_t>(value.internal_to_input.size()), 0 };
    result->effective_bond_orders = { value.effective_bond_orders.data(), static_cast<uint32_t>(value.effective_bond_orders.size()), 0 };
    result->bond_displays = { value.bond_displays.data(), static_cast<uint32_t>(value.bond_displays.size()), 0 };
    result->atom_stereo = { value.atom_stereo.data(), static_cast<uint32_t>(value.atom_stereo.size()), 0 };
    result->clean_pose = value.clean_pose;
}

coordgen_error_t generate(const coordgen_input_t* input, Generation& output, bool collect_probe) {
    const coordgen_error_t input_error = validateInput(input);
    if (input_error != COORDGEN_OK) return input_error;

    PendingMolecule molecule;
    std::vector<sketcherMinimizerAtom*> atoms;
    atoms.reserve(input->atoms.len);
    for (uint32_t index = 0; index < input->atoms.len; ++index) {
        const coordgen_atom_input_t& in = input->atoms.ptr[index];
        sketcherMinimizerAtom* atom = molecule->addNewAtom();
        atom->m_chmN = static_cast<int>(index);
        atom->atomicNumber = static_cast<int>(in.atomic_number);
        atom->charge = in.formal_charge;
        atom->fixed = (in.flags & COORDGEN_ATOM_FIXED) != 0;
        atom->constrained = (in.flags & COORDGEN_ATOM_CONSTRAINED) != 0;
        atom->hidden = (in.flags & COORDGEN_ATOM_HIDDEN) != 0;
        if ((in.flags & COORDGEN_ATOM_HAS_TEMPLATE_COORDINATES) != 0) {
            atom->templateCoordinates = sketcherMinimizerPointF(in.template_coordinates.x, in.template_coordinates.y);
        }
        if ((in.flags & COORDGEN_ATOM_HAS_COORDINATES_3D) != 0) {
            atom->m_x3D = in.coordinates_3d.x;
            atom->m_y3D = in.coordinates_3d.y;
            atom->m_z3D = in.coordinates_3d.z;
        }
        atoms.push_back(atom);
    }
    for (uint32_t index = 0; index < input->atoms.len; ++index) {
        const coordgen_atom_input_t& in = input->atoms.ptr[index];
        if (in.stereo == COORDGEN_ATOM_STEREO_UNSPECIFIED) continue;
        sketcherMinimizerAtomChiralityInfo info;
        info.lookingFrom = in.stereo_looking_from == COORDGEN_INVALID_INDEX ? nullptr : atoms[in.stereo_looking_from];
        info.atom1 = in.stereo_atom_a == COORDGEN_INVALID_INDEX ? nullptr : atoms[in.stereo_atom_a];
        info.atom2 = in.stereo_atom_b == COORDGEN_INVALID_INDEX ? nullptr : atoms[in.stereo_atom_b];
        info.direction = in.stereo == COORDGEN_ATOM_STEREO_COUNTER_CLOCKWISE
                             ? sketcherMinimizerAtomChiralityInfo::counterClockwise
                             : sketcherMinimizerAtomChiralityInfo::clockwise;
        atoms[index]->setStereoChemistry(info);
    }
    std::vector<sketcherMinimizerBond*> bonds;
    bonds.reserve(input->bonds.len);
    for (uint32_t index = 0; index < input->bonds.len; ++index) {
        const coordgen_bond_input_t& in = input->bonds.ptr[index];
        sketcherMinimizerBond* bond = molecule->addNewBond(atoms[in.start], atoms[in.end]);
        bond->m_chmN = static_cast<int>(index);
        bond->bondOrder = static_cast<int>(in.order);
        bond->skip = (in.flags & COORDGEN_BOND_SKIP) != 0;
        bond->crossingBondPenaltyMultiplier = in.crossing_penalty_multiplier;
        bond->hasStereochemistryDisplay = in.display != COORDGEN_BOND_DISPLAY_NONE;
        bond->isWedge = in.display == COORDGEN_BOND_DISPLAY_SOLID_FORWARD || in.display == COORDGEN_BOND_DISPLAY_SOLID_REVERSE;
        bond->isReversed = in.display == COORDGEN_BOND_DISPLAY_SOLID_REVERSE || in.display == COORDGEN_BOND_DISPLAY_HASHED_REVERSE;
        if (in.stereo != COORDGEN_BOND_STEREO_UNSPECIFIED) {
            sketcherMinimizerBondStereoInfo stereo;
            stereo.atom1 = in.stereo_atom_a == COORDGEN_INVALID_INDEX ? nullptr : atoms[in.stereo_atom_a];
            stereo.atom2 = in.stereo_atom_b == COORDGEN_INVALID_INDEX ? nullptr : atoms[in.stereo_atom_b];
            stereo.stereo = (in.stereo == COORDGEN_BOND_STEREO_CIS || in.stereo == COORDGEN_BOND_STEREO_Z)
                                ? sketcherMinimizerBondStereoInfo::cis
                                : sketcherMinimizerBondStereoInfo::trans;
            bond->setStereoChemistry(stereo);
        }
        bonds.push_back(bond);
    }

    sketcherMinimizer minimizer(input->options.precision);
    minimizer.setScoreResidueInteractions(input->options.score_residue_interactions != 0);
    minimizer.setTreatNonterminalBondsToMetalAsZOBs(input->options.treat_nonterminal_bonds_to_metal_as_zero_order != 0);
    minimizer.setEvenAngles(input->options.even_angles != 0);
    minimizer.setSkipMinimization(input->options.skip_minimization != 0);
    minimizer.setForceOpenMacrocycles(input->options.force_open_macrocycles != 0);
    /* build_from_fragments is rejected in validateInput() when nonzero: see
     * the comment there for why there is no faithful call to make here. */
    minimizer.initialize(molecule.release());
    if (input->options.constrain_all_atoms != 0) minimizer.constrainAllAtoms();
    OracleCapture capture;
    CaptureScope capture_scope(collect_probe ? &capture : nullptr);
    output.clean_pose = minimizer.runGenerateCoordinates() ? 1U : 0U;

    output.coordinates.assign(input->atoms.len, { 0.0f, 0.0f });
    output.input_to_internal.assign(input->atoms.len, COORDGEN_INVALID_INDEX);
    output.internal_to_input.reserve(minimizer.getAtoms().size());
    output.atom_stereo.assign(input->atoms.len, COORDGEN_ATOM_STEREO_UNSPECIFIED);
    for (uint32_t index = 0; index < input->atoms.len; ++index) {
        if (!atoms[index]->hidden) {
            const sketcherMinimizerPointF& point = atoms[index]->getCoordinates();
            output.coordinates[index] = { point.x(), point.y() };
            output.atom_stereo[index] = atomStereo(atoms[index]);
        }
    }
    for (uint32_t index = 0; index < minimizer.getAtoms().size(); ++index) {
        const uint32_t original = atomIndex(minimizer.getAtoms()[index]);
        output.internal_to_input.push_back(original);
        if (original != COORDGEN_INVALID_INDEX) output.input_to_internal[original] = index;
    }
    output.effective_bond_orders.reserve(bonds.size());
    output.bond_displays.reserve(bonds.size());
    for (const sketcherMinimizerBond* bond : bonds) {
        output.effective_bond_orders.push_back(static_cast<uint32_t>(bond->bondOrder));
        output.bond_displays.push_back(bondDisplay(bond));
    }
    if (!collect_probe) return COORDGEN_OK;

    std::unordered_map<const sketcherMinimizerRing*, uint32_t> ring_indices;
    std::unordered_map<const sketcherMinimizerFragment*, uint32_t> fragment_indices;
    std::unordered_map<const sketcherMinimizerMolecule*, uint32_t> component_indices;
    std::unordered_map<const sketcherMinimizerMolecule*, const ComponentCoordinates*> before_components;
    std::unordered_map<const sketcherMinimizerMolecule*, const ComponentCoordinates*> after_components;
    for (const ComponentCoordinates& component : capture.before) before_components.emplace(component.component, &component);
    for (const ComponentCoordinates& component : capture.after) after_components.emplace(component.component, &component);
    for (const sketcherMinimizerMolecule* component : minimizer.getMolecules()) {
        const uint32_t component_index = static_cast<uint32_t>(output.components.size());
        component_indices.emplace(component, component_index);
        coordgen_probe_component_t record{};
        record.atom_start = static_cast<uint32_t>(output.component_atoms.size());
        record.atom_count = static_cast<uint32_t>(component->_atoms.size());
        record.transform_status = COORDGEN_PROBE_TRANSFORM_UNOBSERVED;
        const auto before = before_components.find(component);
        const auto after = after_components.find(component);
        setObservedTransform(before == before_components.end() ? nullptr : before->second,
                             after == after_components.end() ? nullptr : after->second,
                             record);
        for (const sketcherMinimizerAtom* atom : component->_atoms) output.component_atoms.push_back(atomIndex(atom));
        output.components.push_back(record);
        std::vector<int> scores;
        /* morganScores indexes its working arrays by each atom's
         * _generalUseN and only documents the requirement in a comment
         * ("assuming that _generalUseN is set as the index of each atom").
         * Fragment building leaves stale values behind, so reading morgan
         * ranks after generation without re-establishing the precondition
         * writes outside the arrays and corrupts the heap. Upstream's own
         * caller, canonicalOrdering(), sets it the same way immediately
         * before the call. */
        for (std::size_t position = 0; position < component->_atoms.size(); ++position) {
            component->_atoms[position]->_generalUseN = static_cast<int>(position);
        }
        const int iterations = sketcherMinimizer::morganScores(component->_atoms, component->_bonds, scores);
        (void)iterations;
        for (int score : scores) output.morgan_ranks.push_back(score < 0 ? 0U : static_cast<uint32_t>(score));
        for (const sketcherMinimizerRing* ring : component->_rings) {
            ring_indices.emplace(ring, static_cast<uint32_t>(output.rings.size()));
            coordgen_probe_ring_t record_ring{ static_cast<uint32_t>(output.ring_atoms.size()), static_cast<uint32_t>(ring->_atoms.size()) };
            for (const sketcherMinimizerAtom* atom : ring->_atoms) output.ring_atoms.push_back(atomIndex(atom));
            output.rings.push_back(record_ring);
        }
    }
    for (uint32_t index = 0; index < minimizer._fragments.size(); ++index) {
        fragment_indices.emplace(minimizer._fragments[index], index);
    }
    std::unordered_map<const sketcherMinimizerFragment*, const TemplateCapture*> template_captures;
    for (const TemplateCapture& template_capture : capture.templates) {
        template_captures[template_capture.fragment] = &template_capture;
    }
    for (sketcherMinimizerFragment* fragment : minimizer._fragments) {
        coordgen_probe_fragment_t record{};
        const auto parent = fragment_indices.find(fragment->getParent());
        record.parent = parent == fragment_indices.end() ? COORDGEN_INVALID_INDEX : parent->second;
        const std::vector<sketcherMinimizerAtom*> fragment_atoms = fragment->getAtoms();
        record.component = fragment_atoms.empty() ? COORDGEN_INVALID_INDEX : component_indices.at(fragment_atoms.front()->getMolecule());
        record.atom_start = static_cast<uint32_t>(output.fragment_atoms.size());
        record.atom_count = static_cast<uint32_t>(fragment_atoms.size());
        for (const sketcherMinimizerAtom* atom : fragment_atoms) output.fragment_atoms.push_back(atomIndex(atom));
        const std::vector<sketcherMinimizerRing*> fragment_rings = fragment->getRings();
        record.ring_start = static_cast<uint32_t>(output.fragment_rings.size());
        record.ring_count = static_cast<uint32_t>(fragment_rings.size());
        for (const sketcherMinimizerRing* ring : fragment_rings) output.fragment_rings.push_back(ring_indices.at(ring));
        record.dof_start = static_cast<uint32_t>(output.dofs.size());
        for (CoordgenFragmentDOF* dof : fragment->getDofs()) {
            coordgen_probe_dof_t probe{};
            probe.id = static_cast<uint32_t>(output.dofs.size());
            probe.kind = dofKind(dof);
            probe.fragment = static_cast<uint32_t>(output.fragments.size());
            probe.current_state = dof->getCurrentState();
            probe.optimal_state = dof->m_optimalState;
            probe.state_count = static_cast<uint32_t>(dof->numberOfStates());
            probe.tier = static_cast<uint32_t>(dof->tier());
            probe.atom_a = COORDGEN_INVALID_INDEX;
            probe.atom_b = COORDGEN_INVALID_INDEX;
            probe.ring = COORDGEN_INVALID_INDEX;
            probe.current_penalty = dof->getCurrentPenalty();
            output.dofs.push_back(probe);
        }
        record.dof_count = static_cast<uint32_t>(output.dofs.size()) - record.dof_start;
        record.flags = (fragment->fixed ? COORDGEN_PROBE_FRAGMENT_FIXED : 0U) |
                       (fragment->isTemplated ? COORDGEN_PROBE_FRAGMENT_TEMPLATED : 0U) |
                       (fragment->constrained ? COORDGEN_PROBE_FRAGMENT_CONSTRAINED : 0U) |
                       (fragment->constrainedFlip ? COORDGEN_PROBE_FRAGMENT_CONSTRAINED_FLIP : 0U) |
                       (fragment->isChain ? COORDGEN_PROBE_FRAGMENT_CHAIN : 0U);
        const auto template_capture = template_captures.find(fragment);
        record.template_match = template_capture == template_captures.end()
                                    ? COORDGEN_INVALID_INDEX
                                    : template_capture->second->template_index;
        record.template_mapping_start = static_cast<uint32_t>(output.template_mapping.size());
        if (template_capture != template_captures.end()) {
            output.template_mapping.insert(output.template_mapping.end(),
                                           template_capture->second->mapping.begin(),
                                           template_capture->second->mapping.end());
        }
        record.template_mapping_count = static_cast<uint32_t>(output.template_mapping.size()) - record.template_mapping_start;
        output.fragments.push_back(record);
    }
    return COORDGEN_OK;
}

void fillProbe(Generation& value, coordgen_probe_result_t* result) {
    result->input_to_internal = { value.input_to_internal.data(), static_cast<uint32_t>(value.input_to_internal.size()), 0 };
    result->internal_to_input = { value.internal_to_input.data(), static_cast<uint32_t>(value.internal_to_input.size()), 0 };
    result->morgan_ranks = { value.morgan_ranks.data(), static_cast<uint32_t>(value.morgan_ranks.size()), 0 };
    result->ring_atoms = { value.ring_atoms.data(), static_cast<uint32_t>(value.ring_atoms.size()), 0 };
    result->fragment_atoms = { value.fragment_atoms.data(), static_cast<uint32_t>(value.fragment_atoms.size()), 0 };
    result->fragment_rings = { value.fragment_rings.data(), static_cast<uint32_t>(value.fragment_rings.size()), 0 };
    result->component_atoms = { value.component_atoms.data(), static_cast<uint32_t>(value.component_atoms.size()), 0 };
    result->template_mapping = const_cast<coordgen_probe_template_mapping_t*>(value.template_mapping.data());
    result->template_mapping_count = static_cast<uint32_t>(value.template_mapping.size());
    result->rings = const_cast<coordgen_probe_ring_t*>(value.rings.data());
    result->ring_count = static_cast<uint32_t>(value.rings.size());
    result->fragments = const_cast<coordgen_probe_fragment_t*>(value.fragments.data());
    result->fragment_count = static_cast<uint32_t>(value.fragments.size());
    result->dofs = const_cast<coordgen_probe_dof_t*>(value.dofs.data());
    result->dof_count = static_cast<uint32_t>(value.dofs.size());
    result->components = const_cast<coordgen_probe_component_t*>(value.components.data());
    result->component_count = static_cast<uint32_t>(value.components.size());
    result->clean_pose = value.clean_pose;
}

} // namespace

extern "C" coordgen_error_t coordgen_generate(const coordgen_input_t* input, coordgen_result_t* result) {
    if (result == nullptr) return COORDGEN_ERROR_INTERNAL;
    std::memset(result, 0, sizeof(*result));
    try {
        auto* owner = new Generation();
        const coordgen_error_t error = generate(input, *owner, false);
        if (error != COORDGEN_OK) { delete owner; return error; }
        fillStable(*owner, result);
        result->owner = owner;
        return COORDGEN_OK;
    } catch (const std::bad_alloc&) {
        return COORDGEN_ERROR_OUT_OF_MEMORY;
    } catch (...) {
        return COORDGEN_ERROR_INTERNAL;
    }
}

extern "C" void coordgen_result_free(coordgen_result_t* result) {
    if (result == nullptr) return;
    delete static_cast<Generation*>(result->owner);
    std::memset(result, 0, sizeof(*result));
}

extern "C" coordgen_error_t coordgen_probe_generate(const coordgen_input_t* input, coordgen_probe_result_t* result) {
    if (result == nullptr) return COORDGEN_ERROR_INTERNAL;
    std::memset(result, 0, sizeof(*result));
    try {
        auto* owner = new Generation();
        const coordgen_error_t error = generate(input, *owner, true);
        if (error != COORDGEN_OK) { delete owner; return error; }
        fillProbe(*owner, result);
        result->owner = owner;
        return COORDGEN_OK;
    } catch (const std::bad_alloc&) {
        return COORDGEN_ERROR_OUT_OF_MEMORY;
    } catch (...) {
        return COORDGEN_ERROR_INTERNAL;
    }
}

extern "C" void coordgen_probe_result_free(coordgen_probe_result_t* result) {
    if (result == nullptr) return;
    delete static_cast<Generation*>(result->owner);
    std::memset(result, 0, sizeof(*result));
}
