#include "moleculesystem.h"
#include "generator.h"

#include <iostream>

using namespace std;

int main()
{
    int nCells = 20;
    double b = 5.620;
    vector<Molecule*> molecules = Generator::generateFcc(nCells, b, AtomType::argon());
    Generator::boltzmannDistributeVelocities(molecules);
    MoleculeSystem system;
    system.addMolecules(molecules);
    system.setBoundaries(0, b*nCells);
    system.setupCells(b*1);

    system.simulate();

    return 0;
}

