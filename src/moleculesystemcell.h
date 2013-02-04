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

    void setIndices(const irowvec& indices);
    const irowvec& indices() const;

    void updateForces();
    void clearMolecules();
private:
    mat geometry;
    vector<MoleculeSystemCell*> m_neighborCells;
    vector<rowvec> m_neighborOffsets;

    mat m_boundaries;
    int m_nDimensions;
    int pow3nDimensions;
    mat cellShiftVectors;

    vector<Molecule*> m_molecules;
    vector<Atom*> m_atoms;

    irowvec m_indices;

    MoleculeSystem* moleculeSystem;
};

#endif // MOLECULESYSTEMCELL_H
