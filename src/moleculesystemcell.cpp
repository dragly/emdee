#include <src/atom.h>

#include <src/moleculesystemcell.h>
#include <src/moleculesystem.h>
#include <src/force/twoparticleforce.h>
#include <src/force/singleparticleforce.h>

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


void MoleculeSystemCell::updateForces()
{
    // Single particle forces
    for(SingleParticleForce* singleParticleForce : m_moleculeSystem->singleParticleForces()) {
        for(Atom* atom : m_atoms) {
            singleParticleForce->apply(atom);
        }
    }

    //    cout << "I have " << m_neighborCells.size() << " neighbors" << endl;
    // Loop over neighbors and their atoms
    for(TwoParticleForce* twoParticleForce : m_moleculeSystem->twoParticleForces()) {
        for(uint iNeighbor = 0; iNeighbor < m_neighborCells.size(); iNeighbor++) {
            MoleculeSystemCell* neighbor = m_neighborCells[iNeighbor];
            if(m_isOnProcessorEdge || neighbor->isOnProcessorEdge()) {
                twoParticleForce->setNewtonsThirdLawEnabled(false);
            } else {
                twoParticleForce->setNewtonsThirdLawEnabled(true);
            }
            const Vector3& neighborOffset = m_neighborOffsets[iNeighbor];
            if(twoParticleForce->isNewtonsThirdLawEnabled()) {
                const irowvec& direction = m_neighborDirections[iNeighbor];
                if(
                        !( // if not one of ..
                           ((direction(0) >= 0 && direction(1) >= 0) && !(direction(0) == 0 && direction(1) == 0 && direction(2) == -1)) // 2x2 in upper right (except right down)
                           || (direction(0) == 1 && direction(1) == -1)) // 1x1 lower right
                        ) // then continue ...
                    //                neighbor->hasAlreadyCalculatedForcesBetweenSelfAndNeighbors()
                    //            )
                {
                    continue;
                }
            }
            const vector<Atom*>& neighborAtoms = neighbor->atoms();
            for(Atom* atom1 : m_atoms) {
                for(Atom* atom2 : neighborAtoms) {
                    if(atom1->isPositionFixed() && atom2->isPositionFixed()) {
                        continue;
                    }
                    twoParticleForce->calculateAndApplyForce(atom1, atom2, neighborOffset);
                }
            }
        }

        // Loop over own atoms
        twoParticleForce->setNewtonsThirdLawEnabled(true);
        for(uint iAtom = 0; iAtom < m_atoms.size(); iAtom++) {
            Atom* atom1 = m_atoms[iAtom];
            int jAtomStart = 0;
            if(twoParticleForce->isNewtonsThirdLawEnabled()) {
                jAtomStart = iAtom + 1;
            }
            for(uint jAtom = jAtomStart; jAtom < m_atoms.size(); jAtom++) {
                if(!twoParticleForce->isNewtonsThirdLawEnabled()) {
                    if(iAtom == jAtom) {
                        continue;
                    }
                }
                Atom* atom2 = m_atoms[jAtom];
                if(atom1->isPositionFixed() && atom2->isPositionFixed()) {
                    continue;
                }
                twoParticleForce->calculateAndApplyForce(atom1, atom2);
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
