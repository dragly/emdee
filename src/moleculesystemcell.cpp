#include "molecule.h"
#include "atom.h"

#include "moleculesystemcell.h"

MoleculeSystemCell::MoleculeSystemCell() :
    m_nDimensions(3),
    pow3nDimensions(pow(3, m_nDimensions)),
    m_indices(zeros<irowvec>(m_nDimensions))
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

void MoleculeSystemCell::addNeighbor(MoleculeSystemCell *cell)
{
    m_neighborCells.push_back(cell);
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

void MoleculeSystemCell::addMolecule(Molecule *molecule) {
    m_molecules.push_back(molecule);
    for(Atom* atom : molecule->atoms()) {
        m_atoms.push_back(atom);
    }
}

void MoleculeSystemCell::setIndices(const irowvec& indices)
{
    m_indices = indices;
}

const irowvec &MoleculeSystemCell::indices() const
{
    return m_indices;
}
