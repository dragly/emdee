#ifndef MOLECULESYSTEMCELL_H
#define MOLECULESYSTEMCELL_H

class Molecule;
class Atom;
class MoleculeSystem;

#include <armadillo>
#include <vector>
#include <iostream>

using namespace arma;
using namespace std;

class MoleculeSystemCell
{
public:
    MoleculeSystemCell(MoleculeSystem *parent);

    friend ostream& operator<<(ostream& os, const MoleculeSystemCell& dt);
    friend ostream& operator<<(ostream& os, MoleculeSystemCell* dt);

    void setBoundaries(mat boundaries);

    void addNeighbor(MoleculeSystemCell *cell, const rowvec& offset);

    const mat &boundaries() const;
    void addMolecule(Molecule *molecule);
    const vector<Atom*>& atoms();
    const vector<Molecule*>& molecules();

    void setIndices(const irowvec& indices);
    const irowvec& indices() const;

    void updateForces();
    void clearMolecules();
//    void calculateForceBetweenAtoms(Atom *atom1, Atom *atom2, const rowvec &neighborOffset, rowvec &rVec, rowvec &force, double& sigma, double& eps);
//    void addAlreadyCalculatedNeighbor(MoleculeSystemCell* neighbor, const rowvec &offset);
    void clearAlreadyCalculatedNeighbors();
    bool hasAlreadyCalculatedForcesBetweenSelfAndNeighbors();
    void setID(int id);
    int id();
private:
    mat geometry;
    vector<MoleculeSystemCell*> m_neighborCells;
//    vector<MoleculeSystemCell*> m_alreadyCalculatedNeighbors;
    vector<rowvec> m_neighborOffsets;
//    vector<rowvec> m_neighborWithAlreadyCalculatedForcesOffsets;

    MoleculeSystem* moleculeSystem;
    bool m_hasAlreadyCalculatedForcesBetweenSelfAndNeighbors;

    mat m_boundaries;
    int m_nDimensions;
    int pow3nDimensions;
    mat cellShiftVectors;

    vector<Molecule*> m_molecules;
    vector<Atom*> m_atoms;

    irowvec m_indices;

    int m_id;
    rowvec force;
};

#endif // MOLECULESYSTEMCELL_H
