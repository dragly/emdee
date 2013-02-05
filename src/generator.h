#ifndef GENERATOR_H
#define GENERATOR_H

class Molecule;
class Atom;

#include "atomtype.h"

#include <vector>
#include <armadillo>
#include <libconfig.h++>

using namespace std;
using namespace arma;
using namespace libconfig;

class Generator
{
public:
    Generator();
    void loadConfiguration(Config* config);
    vector<Molecule*> generateFcc(double b, int nCells, AtomType atomType);
    void boltzmannDistributeVelocities(vector<Molecule *> molecules);

    void setUnitLength(double unitLength);

    const mat& lastBoundaries() const;

    void setNCells(int nCells);
    void setB(double b);
private:
    double m_unitLength;

    mat m_lastBoundaries;
};

#endif // GENERATOR_H
