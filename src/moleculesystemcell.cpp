#include "molecule.h"
#include "atom.h"

#include "moleculesystemcell.h"
#include "moleculesystem.h"

MoleculeSystemCell::MoleculeSystemCell(MoleculeSystem *parent) :
    m_nDimensions(3),
    pow3nDimensions(pow(3, m_nDimensions)),
    m_indices(zeros<irowvec>(m_nDimensions)),
    moleculeSystem(parent)
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

void MoleculeSystemCell::addMolecule(Molecule *molecule) {
    m_molecules.push_back(molecule);
    for(Atom* atom : molecule->atoms()) {
        m_atoms.push_back(atom);
    }
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
    for(Molecule* molecule : m_molecules) {
        molecule->clearForces();
    }
    //    double kB = boltzmannConstant;
    double eps = moleculeSystem->potentialConstant();
    double sigma = 4.5;
    Atom* atom1;
    Atom* atom2;
    rowvec rVec;
    rowvec otherPosition;
//    rowvec shortestVec;
    rowvec cellShiftVector;
//    rowvec cellShiftVectorInUse;
//    double shortestVecSquaredLength;
    double r;
    double sigmar;
    double factor;
    rowvec force;
//    cout << "Updating forces..." << endl;
//    int nCalculations = 0;
//    cout << "m_molecules.size()" << endl;
//    cout << m_molecules.size() << endl;
//    cout << "m_neighborCells.size()" << endl;
//    cout << m_neighborCells.size() << endl;
    for(uint iAtom = 0; iAtom < m_atoms.size(); iAtom++) {
        // Loop over all neighbors (including ourselves)
        for(uint iNeighbor = 0; iNeighbor < m_neighborCells.size(); iNeighbor++) {
            MoleculeSystemCell* neighbor = m_neighborCells.at(iNeighbor);
            rowvec neighborOffset = m_neighborOffsets.at(iNeighbor);
//            cout << "neighbor->atoms().size()" << endl;
//            cout << neighbor->atoms().size() << endl;
            for(uint jAtom = 0; jAtom < neighbor->atoms().size(); jAtom++) {
                atom1 = m_atoms.at(iAtom);
                atom2 = neighbor->atoms().at(jAtom);
                if(atom1 == atom2 && neighbor == this && neighborOffset.max() == 0) {
                    continue;
                }
                cellShiftVector = neighborOffset;
                otherPosition = atom2->absolutePosition() + cellShiftVector;
                rVec = otherPosition - atom1->absolutePosition();
//                rVec = atom2->absolutePosition() + neighborOffset - atom1->absolutePosition();
//                cout << atom2->absolutePosition() << endl;
//                cout << atom1->absolutePosition() << endl;
//                cout << neighborOffset << endl;
//                cout << rVec << endl;
                //            // Minimum image convention
//                shortestVecSquaredLength = INFINITY;
//                for(uint iShiftVec = 0; iShiftVec < moleculeSystem->cellShiftVectors().n_rows; iShiftVec++) {
//                    cellShiftVector = moleculeSystem->cellShiftVectors().row(iShiftVec);
//                    otherPosition = atom2->absolutePosition() + cellShiftVector;
//                    rVec = otherPosition - atom1->absolutePosition();
//                    double rVecSquaredLength = dot(rVec, rVec);
//                    if(rVecSquaredLength < shortestVecSquaredLength) {
//                        shortestVecSquaredLength = rVecSquaredLength;
//                        shortestVec = rVec;
//                        cellShiftVectorInUse = cellShiftVector;
//                    }
//                }
                // Check distances to the nearby cells
                r = norm(rVec, 2);
                sigmar = sigma/r;
                // TODO Verify force term
                factor = - ((24 * eps) / (r*r)) * (2 * pow((sigmar), 12) - pow((sigmar), 6));

                force = factor * rVec;
                atom1->addForce(force);
                //                atom2->addForce(-force);
//                nCalculations++;
            }
        }
    }
//    cout << "nCalculations" << endl;
    //    cout << nCalculations << endl;
}

void MoleculeSystemCell::clearMolecules()
{
    m_molecules.clear();
    m_atoms.clear();
}
