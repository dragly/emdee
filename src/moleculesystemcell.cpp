#include <src/atom.h>

#include <src/moleculesystemcell.h>
#include <src/moleculesystem.h>
#include <src/force/twoparticleforce.h>

MoleculeSystemCell::MoleculeSystemCell(MoleculeSystem *parent) :
    m_nDimensions(3),
    pow3nDimensions(pow(3, m_nDimensions)),
    m_indices(zeros<irowvec>(m_nDimensions)),
    moleculeSystem(parent),
    m_hasAlreadyCalculatedForcesBetweenSelfAndNeighbors(false),
    m_id(0),
    force(blankForce)
{
    cellShiftVectors = zeros(pow3nDimensions, m_nDimensions);
}

void MoleculeSystemCell::setBoundaries(mat boundaries)
{
    m_boundaries = boundaries;
    irowvec counters = zeros<irowvec>(m_nDimensions);
    for(int i = 0; i < pow3nDimensions; i++) {
        irowvec direction = counters - ones<irowvec>(m_nDimensions);
        rowvec lengths = m_boundaries.row(1) - m_boundaries.row(0);
        rowvec directionVec = conv_to<rowvec>::from(direction);
        cellShiftVectors.row(i) = lengths % directionVec;
        counters(0) += 1;
        for(uint idim = 1; idim < counters.size(); idim++) {
            if(counters(idim - 1) > 2) {
                counters(idim - 1) = 0;
                counters(idim) += 1;
            }
        }
    }
    //    cout << cellShiftVectors << endl;
}

void MoleculeSystemCell::addNeighbor(MoleculeSystemCell *cell, const rowvec &offset)
{
    m_neighborCells.push_back(cell);
    m_neighborOffsets.push_back(offset);
}

const mat &MoleculeSystemCell::boundaries() const
{
    return m_boundaries;
}

ostream& operator<<(ostream& os, MoleculeSystemCell* dt)
{
    os << dt->m_boundaries << endl;
    os << dt->m_indices << endl;
    return os;
}

ostream& operator<<(ostream& os, const MoleculeSystemCell& dt)
{
    os << dt.m_boundaries << endl;
    os << dt.m_indices << endl;
    return os;
}

void MoleculeSystemCell::addAtom(Atom *atom) {
    m_atoms.push_back(atom);
    atom->setCellID(m_id);
}

const vector<Atom *> &MoleculeSystemCell::atoms()
{
    return m_atoms;
}

void MoleculeSystemCell::setIndices(const irowvec& indices)
{
    m_indices = indices;
}

const irowvec &MoleculeSystemCell::indices() const
{
    return m_indices;
}


void MoleculeSystemCell::updateForces()
{
    TwoParticleForce* interatomicForce = moleculeSystem->interatomicForce();

    // Loop over neighbors and their atoms
    for(uint iNeighbor = 0; iNeighbor < m_neighborCells.size(); iNeighbor++) {
        MoleculeSystemCell* neighbor = m_neighborCells[iNeighbor];
        const rowvec& neighborOffset = m_neighborOffsets[iNeighbor];
        if(neighbor->hasAlreadyCalculatedForcesBetweenSelfAndNeighbors()) {
            continue;
        }
        const vector<Atom*>& neighborAtoms = neighbor->atoms();
        for(Atom* atom1 : m_atoms) {
            for(Atom* atom2 : neighborAtoms) {
                interatomicForce->calculateAndApplyForce(atom1, atom2, neighborOffset);
            }
        }
        m_hasAlreadyCalculatedForcesBetweenSelfAndNeighbors = true;
    }

    // Loop over own atoms
    for(uint iAtom = 0; iAtom < m_atoms.size(); iAtom++) {
        Atom* atom1 = m_atoms[iAtom];
        for(uint jAtom = iAtom + 1; jAtom < m_atoms.size(); jAtom++) {
            Atom* atom2 = m_atoms[jAtom];
            interatomicForce->calculateAndApplyForce(atom1, atom2);
        }
    }
}

void MoleculeSystemCell::clearAtoms()
{
    m_atoms.clear();
}

void MoleculeSystemCell::clearAlreadyCalculatedNeighbors()
{
    m_hasAlreadyCalculatedForcesBetweenSelfAndNeighbors = false;
}

void MoleculeSystemCell::setID(int id)
{
    m_id = id;
}

int MoleculeSystemCell::id()
{
    return m_id;
}
