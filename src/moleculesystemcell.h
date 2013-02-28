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
    void clearAlreadyCalculatedNeighbors();
    void setID(int id);
    int id();

    bool hasAlreadyCalculatedForcesBetweenSelfAndNeighbors()
    {
        return m_hasAlreadyCalculatedForcesBetweenSelfAndNeighbors;
    }
protected:
    mat geometry;
    vector<MoleculeSystemCell*> m_neighborCells;
    vector<rowvec> m_neighborOffsets;

    int m_nDimensions;
    int pow3nDimensions;
    irowvec m_indices;
    MoleculeSystem* moleculeSystem;
    bool m_hasAlreadyCalculatedForcesBetweenSelfAndNeighbors;

    mat m_boundaries;
    mat cellShiftVectors;

    vector<Molecule*> m_molecules;
    vector<Atom*> m_atoms;


    int m_id;
    rowvec& force;
    rowvec blankForce;
};

#endif // MOLECULESYSTEMCELL_H
