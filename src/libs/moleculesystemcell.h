#ifndef MOLECULESYSTEMCELL_H
#define MOLECULESYSTEMCELL_H

class Atom;
class MoleculeSystem;

#include <armadillo>
#include <vector>
#include <iostream>

#include <math/vector3.h>

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

    void updateTwoParticleForceAndNeighborAtoms();
    void clearAtoms();
    void clearAlreadyCalculatedNeighbors();
    void deleteAtomsFromCellAndSystem();
    void setID(int id);
    int id();
    void addAtoms(const vector<Atom *> &atoms);
    void deleteAtoms(int nAtoms);
    const vector<MoleculeSystemCell*> &neighborCells() {
        return m_neighborCells;
    }

    void setLocal(bool localCell);
    bool local() const;
protected:
    mat geometry;
    vector<MoleculeSystemCell*> m_neighborCells;
    vector<Vector3> m_neighborOffsets;
    vector<irowvec> m_neighborDirections;

    int m_nDimensions;
    int pow3nDimensions;
    irowvec m_indices;
    MoleculeSystem* m_moleculeSystem;

    mat m_boundaries;
    mat cellShiftVectors;

    vector<Atom*> m_atoms;


    int m_id;
    Vector3& force;
    Vector3 blankForce;
    Vector3 zeroOffset;
    bool m_isOnProcessorEdge;
    bool m_isLocalCell;
};

#endif // MOLECULESYSTEMCELL_H
