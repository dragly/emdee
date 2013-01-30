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
    static vector<Molecule*> generateFcc(int nCells, double b, AtomType atomType);
    static void boltzmannDistributeVelocities(vector<Molecule *> molecules);
};

#endif // GENERATOR_H
