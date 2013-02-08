#include <src/moleculesystem.h>
#include <src/generator.h>
#include <src/interatomicforce.h>
#include <src/integrator/integrator.h>
#include <src/integrator/velocityverletintegrator.h>
#include <src/integrator/eulercromerintegrator.h>

#include <iostream>
#include <libconfig.h++>

using namespace std;
using namespace libconfig;

int main(/*int argc, char** argv*/)
{
    Config config;
    config.readFile("testconfig.cfg");

    double unitLength = config.lookup("units.length");
//    double unitTime = config.lookup("units.time");
    double unitTemperature = config.lookup("units.temperature");

    int nSimulationSteps = config.lookup("simulation.nSimulationSteps");
    double timeStep = config.lookup("integrator.timeStep");
    Generator generator;
    generator.loadConfiguration(&config);

    // Generator specific config
    int nCells = config.lookup("generator.fcc.nCells");
    double b = config.lookup("generator.fcc.b");
    double temperature = config.lookup("system.initialTemperature");
    temperature /= unitTemperature;
    b /= unitLength;
    double potentialConstant = config.lookup("system.potentialConstant");
    potentialConstant /= unitLength;
    vector<Molecule*> molecules = generator.generateFcc(b, nCells, AtomType::argon());
    generator.boltzmannDistributeVelocities(temperature, molecules);
    // Set up force
    InteratomicForce* force = new InteratomicForce();
    force->setPotentialConstant(potentialConstant);
    // Set up molecule system
    MoleculeSystem system;
    system.setInteratomicForce(force);
    Integrator* integrator = new VelocityVerletIntegrator(&system);
    integrator->setTimeStep(timeStep);
    system.setIntegrator(integrator);
    system.loadConfiguration(&config); // TODO remove this
    system.addMolecules(molecules);
    cout << "addded" << endl;
    system.setBoundaries(generator.lastBoundaries());
    cout << "setbounds" << endl;
    system.setupCells(potentialConstant * 3);
    cout << "Setup cells" << endl;
    system.simulate(nSimulationSteps);
    ofstream test;
    test.open("test");
    return 0;
}

