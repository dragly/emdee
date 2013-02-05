#ifndef GENERATOR_H
#define GENERATOR_H

class Molecule;
class Atom;

#include "atomtype.h"

#include <vector>
#include <armadillo>

using namespace std;
using namespace arma;

class Generator
{
public:
    Generator();
    vector<Molecule*> generateFcc(int nCells, double b, AtomType atomType);
    void boltzmannDistributeVelocities(vector<Molecule *> molecules);

    void setUnitLength(double unitLength);

private:
    double m_unitLength;
};

#endif // GENERATOR_H
