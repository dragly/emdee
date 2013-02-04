#ifndef MOLECULESYSTEMCELL_H
#define MOLECULESYSTEMCELL_H

class Molecule;
class Atom;

#include <armadillo>
#include <vector>
#include <iostream>

using namespace arma;
using namespace std;

class MoleculeSystemCell
{
public:
    MoleculeSystemCell();

    friend ostream& operator<<(ostream& os, const MoleculeSystemCell& dt);
    friend ostream& operator<<(ostream& os, MoleculeSystemCell* dt);

    void setBoundaries(mat boundaries);

    void addNeighbor(MoleculeSystemCell *cell);

    const mat &boundaries() const;
    void addMolecule(Molecule *molecule);

    void setIndices(const irowvec& indices);
    const irowvec& indices() const;

    void updateForces();
private:
    mat geometry;
    vector<MoleculeSystemCell*> m_neighborCells;

    mat m_boundaries;
    int m_nDimensions;
    int pow3nDimensions;
    mat cellShiftVectors;

    vector<Molecule*> m_molecules;
    vector<Atom*> m_atoms;

    irowvec m_indices;
};

#endif // MOLECULESYSTEMCELL_H
