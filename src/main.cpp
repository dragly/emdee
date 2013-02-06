#include "moleculesystem.h"
#include "generator.h"

#include <iostream>
#include <libconfig.h++>

using namespace std;
using namespace libconfig;

int main(/*int argc, char** argv*/)
{
    Config config;
    config.readFile("testconfig.cfg");

    double unitLength = config.lookup("units.length");

    int nSimulationSteps = config.lookup("simulation.nSimulationSteps");
    Generator generator;
    generator.loadConfiguration(&config);

    // Generator specific config
    int nCells = config.lookup("generator.fcc.nCells");
    double b = config.lookup("generator.fcc.b");
    double bUnit = b / unitLength;
    double potentialConstant = config.lookup("system.potentialConstant");
    double potentialConstantUnit = potentialConstant / unitLength;
    vector<Molecule*> molecules = generator.generateFcc(bUnit, nCells, AtomType::argon());
    generator.boltzmannDistributeVelocities(molecules);
    MoleculeSystem system;
    system.loadConfiguration(&config);
    system.addMolecules(molecules);
    cout << "addded" << endl;
    system.setBoundaries(generator.lastBoundaries());
    cout << "setbounds" << endl;
    system.setupCells();
    cout << "Setup cells" << endl;
    system.simulate(nSimulationSteps);

    return 0;
}

