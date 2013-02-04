#ifndef MOLECULESYSTEMCELL_H
#define MOLECULESYSTEMCELL_H

#include <armadillo>
#include <vector>
#include <iostream>

using namespace arma;
using namespace std;

class MoleculeSystemCell
{
public:
    MoleculeSystemCell();

    void setBoundaries(mat boundaries);
    friend ostream& operator<<(ostream& os, const MoleculeSystemCell& dt);
    friend ostream& operator<<(ostream& os, MoleculeSystemCell* dt);
private:
    mat geometry;
    vector<MoleculeSystemCell*> neighborCells;

    mat m_boundaries;
    int m_nDimensions;
    int pow3nDimensions;
    mat cellShiftVectors;
};

#endif // MOLECULESYSTEMCELL_H
