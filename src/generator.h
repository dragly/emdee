#ifndef GENERATOR_H
#define GENERATOR_H

class Molecule;
class Atom;

#include <src/atomtype.h>

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

    vector<Molecule*> generateFcc(double sideLength, int nCells, AtomType atomType);
    void boltzmannDistributeVelocities(double temperature, const vector<Molecule *> &molecules);
    void uniformDistributeVelocities(double maxVelocity, vector<Molecule *> molecules);

    const mat& lastBoundaries() const;

    void setNCells(int nCells);
    void setB(double b);
protected:
    double m_unitLength;

    mat m_lastBoundaries;

    int m_nDimensions;
};

#endif // GENERATOR_H
