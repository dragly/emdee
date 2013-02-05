#include "moleculesystem.h"
#include "generator.h"

#include <iostream>

using namespace std;

int main()
{
    double sigma = 3.405;
    int nCells = 6;
    double b = 5.620;
    Generator generator;
    generator.setUnitLength(sigma);
    vector<Molecule*> molecules = generator.generateFcc(nCells, b, AtomType::argon());
    generator.boltzmannDistributeVelocities(molecules);
    MoleculeSystem system;
    system.setUnitLength(sigma);
    system.setPotentialConstant(sigma);
    system.addMolecules(molecules);
    system.setBoundaries(0, b*nCells);
    system.setupCells(sigma*3);

    system.simulate(10000);

    return 0;
}

