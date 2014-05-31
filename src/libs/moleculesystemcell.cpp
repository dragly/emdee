#include <atom.h>

#include <moleculesystemcell.h>
#include <moleculesystem.h>
#include <force/twoparticleforce.h>
#include <force/threeparticleforce.h>
#include <force/singleparticleforce.h>

MoleculeSystemCell::MoleculeSystemCell(MoleculeSystem *parent) :
    m_nDimensions(3),
    pow3nDimensions(pow(3, m_nDimensions)),
    //    m_indices(zeros<iVector3>(m_nDimensions)),
    m_moleculeSystem(parent),
    //    m_hasAlreadyCalculatedForcesBetweenSelfAndNeighbors(false),
    m_id(0),
    force(blankForce),
    m_isOnProcessorEdge(false),
    m_isLocalCell(true)
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
        for(uint idim = 1; idim < 3; idim++) {
            if(counters(idim - 1) > 2) {
                counters(idim - 1) = 0;
                counters(idim) += 1;
            }
        }
    }
    //    cout << cellShiftVectors << endl;
}

void MoleculeSystemCell::addNeighbor(MoleculeSystemCell *cell, const Vector3 &offset, const irowvec& direction)
{
    m_neighborCells.push_back(cell);
    m_neighborOffsets.push_back(offset);
    m_negativeNeighborOffsets.push_back(-offset);
    m_neighborDirections.push_back(direction);
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

void MoleculeSystemCell::setIndices(const irowvec & indices)
{
    m_indices = indices;
}

const irowvec &MoleculeSystemCell::indices() const
{
    return m_indices;
}

void MoleculeSystemCell::updateTwoParticleForceAndNeighborAtoms()
{
    TwoParticleForce* twoParticleForce = m_moleculeSystem->twoParticleForce();
    bool needsNeighborLists = (m_moleculeSystem->threeParticleForce() != NULL);
    // Single particle forces
    for(SingleParticleForce* singleParticleForce : m_moleculeSystem->singleParticleForces()) {
        for(Atom* atom : m_atoms) {
            singleParticleForce->apply(atom);
        }
    }

    // Loop over neighbors and their atoms
    if(twoParticleForce) {
        double cutoffRadiusSquared = twoParticleForce->cutoffRadius()*twoParticleForce->cutoffRadius();
        twoParticleForce->setNewtonsThirdLawEnabled(true);
        uint neighborCellCount = m_neighborCells.size();
        for(uint iNeighbor = 0; iNeighbor < neighborCellCount; iNeighbor++) {
            MoleculeSystemCell* neighborCell = m_neighborCells[iNeighbor];
            // No reason to calculate forces between atoms on two ghost cells.
            // This will also skip internal calculations on ghost cells.
            // This skip may only be performed if there are no need for
            // listing neighbors, as it is when three-particle forces are in use
            if(!m_moleculeSystem->threeParticleForce() && !m_isLocalCell && !neighborCell->local()) {
                continue;
            }
            Vector3 &neighborOffset = (m_neighborOffsets[iNeighbor]);
            double neighborOffsetSquared = Vector3::dot(neighborOffset, neighborOffset);
            Vector3 &negativeNeighborOffset = (m_negativeNeighborOffsets[iNeighbor]);

            // Exception to ID check below is when the neighbor cell is a copy of this,
            // but still a neighbor. This results in a situations where
            // the atoms will have the same IDs, even though we want to
            // calculate their forces
            bool neighborIsCopyOfThis = (this == neighborCell && neighborOffsetSquared > 1e-12);

            const vector<Atom*>& neighborCellAtoms = neighborCell->atoms();
            for(Atom* atom1 : m_atoms) {
                for(Atom* atom2 : neighborCellAtoms) {
                    // Newton's third law is implemented on the next line.
                    // It works so that each atom is responsible for
                    // calculating and applying forces to all atoms with a higher atom ID
                    // This is an alternative to only caring about neighbor cells in a
                    // certain direction.
                    // Note the above exception when the neighbor is the same as this,
                    // but still a neighbor.
                    if(atom1->id() >= atom2->id() && !neighborIsCopyOfThis) {
                        continue;
                    }
                    // Skip if position of either is fixed
                    if(atom1->isPositionFixed() && atom2->isPositionFixed()) {
                        continue;
                    }
                    // Skip if outside cutoff radius
                    double distanceSquared = Vector3::distanceSquared(atom1->position(), atom2->position() + neighborOffset);
                    if(distanceSquared > cutoffRadiusSquared) {
                        continue;
                    }
                    if(needsNeighborLists) {
                        atom1->addNeighborAtom(atom2, &neighborOffset);
                        if(!neighborIsCopyOfThis) {
                            atom2->addNeighborAtom(atom1, &negativeNeighborOffset);
                        }
                    }

                    // No reason to calculate forces between atoms on two ghost cells.
                    // This will also skip internal calculations on ghost cells.
                    // We had to add neighbors above, however
                    if(!m_isLocalCell && !neighborCell->local()) {
                        continue;
                    }
                    // Disable Newton's third law if we are acting on a copy of our own cell
                    if(neighborIsCopyOfThis) {
                        twoParticleForce->setNewtonsThirdLawEnabled(false);
                    } else {
                        twoParticleForce->setNewtonsThirdLawEnabled(true);
                    }
                    twoParticleForce->calculateAndApplyForce(atom1, atom2, neighborOffset);
                }
            }
        }
    }
}

void MoleculeSystemCell::clearAtoms()
{
    m_atoms.clear();
}

void MoleculeSystemCell::addAtoms(const vector<Atom*>& atoms) {
    m_atoms.insert(m_atoms.end(), atoms.begin(), atoms.end());
}

//void MoleculeSystemCell::clearAlreadyCalculatedNeighbors()
//{
//    m_hasAlreadyCalculatedForcesBetweenSelfAndNeighbors = false;
//}

void MoleculeSystemCell::setID(int id)
{
    m_id = id;
}

int MoleculeSystemCell::id() const
{
    return m_id;
}


void MoleculeSystemCell::deleteAtomsFromCellAndSystem()
{
    vector<Atom*> atomsExceptThisCell;
    for(MoleculeSystemCell* cell : m_moleculeSystem->globalCells()) {
        if(cell == this) {
            continue;
        }
        atomsExceptThisCell.insert(atomsExceptThisCell.end(), cell->atoms().begin(), cell->atoms().end());
    }
    m_moleculeSystem->clearAtoms();
    m_moleculeSystem->addAtoms(atomsExceptThisCell);
    for(Atom* atom : m_atoms) {
        delete atom;
    }
    m_atoms.clear();
}

void MoleculeSystemCell::deleteAtoms(int nAtoms) {
    m_atoms.erase(m_atoms.end() - nAtoms, m_atoms.end());
}

void MoleculeSystemCell::setLocal(bool localCell)
{
    m_isLocalCell = localCell;
}

bool MoleculeSystemCell::local() const
{
    return m_isLocalCell;
}
