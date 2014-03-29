#include <atom.h>

#include <moleculesystemcell.h>
#include <moleculesystem.h>
#include <force/twoparticleforce.h>
#include <force/threeparticleforce.h>
#include <force/singleparticleforce.h>

#ifdef MD_USE_GLOG
#include <glog/logging.h>
#else
#include <glogfallback.h>
#endif

MoleculeSystemCell::MoleculeSystemCell(MoleculeSystem *parent) :
    m_nDimensions(3),
    pow3nDimensions(pow(3, m_nDimensions)),
    //    m_indices(zeros<iVector3>(m_nDimensions)),
    m_moleculeSystem(parent),
    //    m_hasAlreadyCalculatedForcesBetweenSelfAndNeighbors(false),
    m_id(0),
    force(blankForce),
    m_isOnProcessorEdge(false)
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

bool MoleculeSystemCell::shouldNewtonsThirdBeEnabled(MoleculeSystemCell* neighbor) {
    (void)neighbor;
    return true;
    //    if(m_isOnProcessorEdge || neighbor->isOnProcessorEdge()) {
    //        return false;
    //    } else {
    //        return true;
    //    }
}

bool MoleculeSystemCell::checkDirection(int neighborID) {
    const irowvec& direction = m_neighborDirections[neighborID];
    return checkDirection(direction);
}

bool MoleculeSystemCell::checkDirection(const irowvec &direction) {
    bool is2x2x3UpperRight = ((direction(0) >= 0 && direction(1) >= 0)
                              && !(direction(0) == 0 && direction(1) == 0 && direction(2) == -1));
    bool is1x1x3LowerRight = (direction(0) == 1 && direction(1) == -1);
    return is2x2x3UpperRight || is1x1x3LowerRight;
}

void MoleculeSystemCell::updateForces()
{
    LOG(INFO) << "Updating forces in cell";
    TwoParticleForce* twoParticleForce = m_moleculeSystem->twoParticleForce();
    ThreeParticleForce* threeParticleForce = m_moleculeSystem->threeParticleForce();
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
        for(uint iNeighbor = 0; iNeighbor < m_neighborCells.size(); iNeighbor++) {
            MoleculeSystemCell* neighborCell = m_neighborCells[iNeighbor];
            const Vector3& neighborOffset = m_neighborOffsets[iNeighbor];

            if(!checkDirection(iNeighbor)) {
                continue;
            }

            const vector<Atom*>& neighborCellAtoms = neighborCell->atoms();
            for(Atom* atom1 : m_atoms) {
                for(Atom* atom2 : neighborCellAtoms) {
                    if(atom1 == atom2) {
                        continue;
                    }
                    if(atom1->isPositionFixed() && atom2->isPositionFixed()) {
                        continue;
                    }
                    double distanceSquared = Vector3::differenceSquared(atom1->position(), atom2->position());
                    if(distanceSquared > cutoffRadiusSquared) {
                        continue;
                    }
                    atom1->addNeighborAtom(atom2, &neighborOffset);
                    atom2->addNeighborAtom(atom1, &neighborOffset);
                    twoParticleForce->calculateAndApplyForce(atom1, atom2, neighborOffset);
                }
            }
        }

        // Loop over own atoms
//        for(uint iAtom = 0; iAtom < m_atoms.size(); iAtom++) {
//            Atom* atom1 = m_atoms[iAtom];
//            int jAtomStart = 0;
//            jAtomStart = iAtom + 1;
//            for(uint jAtom = jAtomStart; jAtom < m_atoms.size(); jAtom++) {
//                Atom* atom2 = m_atoms[jAtom];
//                if(atom1->isPositionFixed() && atom2->isPositionFixed()) {
//                    continue;
//                }
//                double distanceSquared = Vector3::differenceSquared(atom1->position(), atom2->position());
//                if(distanceSquared > cutoffRadiusSquared) {
//                    continue;
//                }
//                atom1->addNeighborAtom(atom2, &zeroOffset);
//                atom2->addNeighborAtom(atom1, &zeroOffset);
//                twoParticleForce->calculateAndApplyForce(atom1, atom2);
//            }
//        }
    }

    if(threeParticleForce) {
        for(uint iAtom = 0; iAtom < m_atoms.size(); iAtom++) {
            Atom* atom1 = m_atoms[iAtom];
            for(uint jAtom = 0; jAtom < atom1->neighborAtoms().size(); jAtom++) {
                std::pair<Atom*,const Vector3*> atom2 = atom1->neighborAtoms()[jAtom];
                for(uint kAtom = jAtom + 1; kAtom < atom1->neighborAtoms().size(); kAtom++) {
                    std::pair<Atom*,const Vector3*> atom3 = atom1->neighborAtoms()[kAtom];
                    const Vector3& offsetVector2 = (*(atom2.second));
                    const Vector3& offsetVector3 = (*(atom3.second));
                    threeParticleForce->calculateAndApplyForce(atom1, atom2.first, atom3.first, offsetVector2, offsetVector3);
                }
            }
        }
    }

//    // TODO Ensure the three-particle forces work properly by creating a test that checks the potential for a run with one big cell and one with 3x3x3 cells

//    // Three particle forces
//    Vector3 zeroVector;
//    zeroVector.zeros();

//    // Two neighbor atoms
//    if(threeParticleForce) {
//        for(uint iNeighbor1 = 0; iNeighbor1 < m_neighborCells.size(); iNeighbor1++) {
//            MoleculeSystemCell* neighbor1 = m_neighborCells[iNeighbor1];
//            const Vector3& neighborOffset1 = m_neighborOffsets[iNeighbor1];
//            const vector<Atom*>& neighborAtoms1 = neighbor1->atoms();

//            threeParticleForce->setNewtonsThirdLawEnabled(shouldNewtonsThirdBeEnabled(neighbor1));
//            if(threeParticleForce->isNewtonsThirdLawEnabled() && !checkDirection(iNeighbor1)) {
//                continue;
//            }
//            int iNeighbor2Start = 0;
//            if(threeParticleForce->isNewtonsThirdLawEnabled()) {
//                iNeighbor2Start = iNeighbor1;
//            }
//            for(uint iNeighbor2 = iNeighbor2Start; iNeighbor2 < m_neighborCells.size(); iNeighbor2++) {
//                MoleculeSystemCell* neighbor2 = m_neighborCells[iNeighbor2];
//                const Vector3& neighborOffset2 = m_neighborOffsets[iNeighbor2];
//                const vector<Atom*>& neighborAtoms2 = neighbor2->atoms();

//                threeParticleForce->setNewtonsThirdLawEnabled(shouldNewtonsThirdBeEnabled(neighbor2));
//                if(threeParticleForce->isNewtonsThirdLawEnabled() && !checkDirection(iNeighbor2)) {
//                    continue;
//                }
//                for(Atom* atom1 : m_atoms) {
//                    for(uint jAtom = 0; jAtom < neighborAtoms1.size(); jAtom++) {
//                        Atom* atom2 = neighborAtoms1[jAtom];
//                        uint kAtomStart = 0;
//                        if(iNeighbor1 == iNeighbor2) {
//                            kAtomStart = jAtom + 1;
//                        }
//                        for(uint kAtom = kAtomStart; kAtom < neighborAtoms2.size(); kAtom++) {
//                            Atom* atom3 = neighborAtoms2[kAtom];
//                            if(atom1->isPositionFixed() && atom2->isPositionFixed() && atom3->isPositionFixed()) {
//                                continue;
//                            }
//                            threeParticleForce->calculateAndApplyForce(atom1, atom2, atom3, neighborOffset1, neighborOffset2);
//                            //                            cout << m_id << endl;
//                            //                            cout << "Two in neighbor!" << endl << endl << endl;
//                        }
//                    }
//                }
//            }
//        }

//        // Two local atoms
//        for(uint iNeighbor1 = 0; iNeighbor1 < m_neighborCells.size(); iNeighbor1++) {
//            MoleculeSystemCell* neighbor1 = m_neighborCells[iNeighbor1];
//            const Vector3& neighborOffset1 = m_neighborOffsets[iNeighbor1];
//            const vector<Atom*>& neighborAtoms1 = neighbor1->atoms();

//            threeParticleForce->setNewtonsThirdLawEnabled(shouldNewtonsThirdBeEnabled(neighbor1));
//            if(threeParticleForce->isNewtonsThirdLawEnabled() && !checkDirection(iNeighbor1)) {
//                continue;
//            }
//            for(uint iAtom = 0; iAtom < m_atoms.size(); iAtom++) {
//                Atom* atom1 = m_atoms[iAtom];
//                int jAtomStart = 0;
//                if(threeParticleForce->isNewtonsThirdLawEnabled()) {
//                    jAtomStart = iAtom + 1;
//                }
//                for(uint jAtom = jAtomStart; jAtom < m_atoms.size(); jAtom++) {
//                    Atom* atom2 = m_atoms[jAtom];
//                    if(!threeParticleForce->isNewtonsThirdLawEnabled()) {
//                        if(atom1 == atom2) {
//                            continue;
//                        }
//                    }
//                    for(Atom* atom3 : neighborAtoms1) {
//                        if(atom1->isPositionFixed() && atom2->isPositionFixed() && atom3->isPositionFixed()) {
//                            continue;
//                        }
//                        threeParticleForce->calculateAndApplyForce(atom1, atom2, atom3, zeroVector, neighborOffset1);
//                        //                        threeParticleForce->calculateAndApplyForce(atom2, atom1, atom3, zeroVector, neighborOffset1);
//                        //                        cout << m_id << endl;
//                        //                        cout << iAtom << " " << jAtom << endl;
//                        //                        cout << "Two in own!" << endl << endl << endl;
//                    }
//                }
//            }
//        }

//        // Loop over own atoms
//        threeParticleForce->setNewtonsThirdLawEnabled(true);
//        for(uint iAtom = 0; iAtom < m_atoms.size(); iAtom++) {
//            Atom* atom1 = m_atoms[iAtom];
//            int jAtomStart = 0;
//            if(threeParticleForce->isNewtonsThirdLawEnabled()) {
//                jAtomStart = iAtom + 1;
//            }
//            for(uint jAtom = jAtomStart; jAtom < m_atoms.size(); jAtom++) {
//                if(!threeParticleForce->isNewtonsThirdLawEnabled()) {
//                    if(iAtom == jAtom) {
//                        continue;
//                    }
//                }
//                Atom* atom2 = m_atoms[jAtom];
//                int kAtomStart = 0;
//                if(threeParticleForce->isNewtonsThirdLawEnabled()) {
//                    kAtomStart = jAtom + 1;
//                }
//                for(uint kAtom = kAtomStart; kAtom < m_atoms.size(); kAtom++) {
//                    if(!threeParticleForce->isNewtonsThirdLawEnabled()) {
//                        if(iAtom == jAtom || jAtom == kAtom || iAtom == kAtom) {
//                            continue;
//                        }
//                    }
//                    Atom* atom3 = m_atoms[kAtom];
//                    if(atom1->isPositionFixed() && atom2->isPositionFixed() && atom3->isPositionFixed()) {
//                        continue;
//                    }
//                    //                    cout << iAtom << " " << jAtom << " " << kAtom << endl;
//                    threeParticleForce->calculateAndApplyForce(atom1, atom2, atom3);
//                    //                    cout << "Force " << atom1->force() << "; " << atom2->force() << "; " << atom3->force() << endl;
//                    //                    cout << m_id << endl;
//                    //                    cout << "All in own!" << endl << endl << endl;
//                }
//            }
//        }
//    }
    //    }
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

int MoleculeSystemCell::id()
{
    return m_id;
}


void MoleculeSystemCell::deleteAtomsFromCellAndSystem()
{
    vector<Atom*> atomsExceptThisCell;
    for(MoleculeSystemCell* cell : m_moleculeSystem->cells()) {
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
