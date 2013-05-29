#ifndef MOLECULESYSTEMCELL_H
#define MOLECULESYSTEMCELL_H

class Atom;
class MoleculeSystem;

#include <armadillo>
#include <vector>
#include <iostream>

#include <src/math/vector3.h>

using namespace arma;
using namespace std;

class MoleculeSystemCell
{
public:
    MoleculeSystemCell(MoleculeSystem *parent);

    friend ostream& operator<<(ostream& os, const MoleculeSystemCell& dt);
    friend ostream& operator<<(ostream& os, MoleculeSystemCell* dt);

    void setBoundaries(mat boundaries);

    void addNeighbor(MoleculeSystemCell *cell, const Vector3& offset, const irowvec &direction);

    const mat &boundaries() const;
    void addAtom(Atom *atom);
    const vector<Atom*>& atoms();

    void setIndices(const irowvec& indices);
    const irowvec& indices() const;

    void updateForces();
    void clearAtoms();
    void clearAlreadyCalculatedNeighbors();
    void deleteAtomsFromCellAndSystem();
    void setID(int id);
    int id();

//    bool hasAlreadyCalculatedForcesBetweenSelfAndNeighbors()
//    {
//        return m_hasAlreadyCalculatedForcesBetweenSelfAndNeighbors;
    //    }
    void addAtoms(const vector<Atom *> &atoms);
    void setOnProcessorEdge(bool enable);
    bool isOnProcessorEdge();
    void deleteAtoms(int nAtoms);
    bool shouldNewtonsThirdBeEnabled(MoleculeSystemCell *neighbor);
    bool checkDirection(int neighborID);
    bool checkDirection(const irowvec &direction);
protected:
    mat geometry;
    vector<MoleculeSystemCell*> m_neighborCells;
    vector<Vector3> m_neighborOffsets;
    vector<irowvec> m_neighborDirections;

    int m_nDimensions;
    int pow3nDimensions;
    irowvec m_indices;
    MoleculeSystem* m_moleculeSystem;
//    bool m_hasAlreadyCalculatedForcesBetweenSelfAndNeighbors;

    mat m_boundaries;
    mat cellShiftVectors;

    vector<Atom*> m_atoms;


    int m_id;
    Vector3& force;
    Vector3 blankForce;
    bool m_isOnProcessorEdge;
};

inline void MoleculeSystemCell::setOnProcessorEdge(bool enable) {
    m_isOnProcessorEdge = enable;
}

inline bool MoleculeSystemCell::isOnProcessorEdge() {
    return m_isOnProcessorEdge;
}

#endif // MOLECULESYSTEMCELL_H
