#include <src/moleculesystem.h>
#include <src/generator.h>
#include <src/interatomicforce.h>
#include <src/integrator/integrator.h>
#include <src/integrator/velocityverletintegrator.h>
#include <src/integrator/eulercromerintegrator.h>
#include <src/filemanager.h>

#include <iostream>
#include <libconfig.h++>

using namespace std;
using namespace libconfig;

int main(/*int argc, char** argv*/)
{
    Config config;
    config.readFile("testconfig.cfg");

    double unitLength = config.lookup("units.length");
    double unitTime = config.lookup("units.time");
    double unitMass = config.lookup("units.mass");
    double unitTemperature = config.lookup("units.temperature");
    double boltzmannConstant = config.lookup("units.boltzmannConstant");
    boltzmannConstant /= (unitLength*unitLength * unitMass / (unitTime*unitTime));

    // File manager config
    string outFileName;
    config.lookupValue("simulation.saveFileName", outFileName);

    int nSimulationSteps = config.lookup("simulation.nSimulationSteps");
    double timeStep = config.lookup("integrator.timeStep");
    timeStep /= unitTime;
    Generator generator;
//    generator.loadConfiguration(&config);

    // Generator specific config
    string velocityDistributionType;
    config.lookupValue("generator.velocity.distribution", velocityDistributionType);
    int nCells = config.lookup("generator.fcc.nCells");
    double b = config.lookup("generator.fcc.b");
    double temperature = config.lookup("system.initialTemperature");
    temperature /= unitTemperature;
    b /= unitLength;

    // Force and potentials
    double potentialConstant = config.lookup("system.potentialConstant");
    potentialConstant /= unitLength;
    double energyConstant = config.lookup("system.energyConstant");
    energyConstant /= (unitMass * unitLength * unitLength / (unitTime * unitTime));

    // Set up grid
    vector<Molecule*> molecules = generator.generateFcc(b, nCells, AtomType::argon());

    // Distribute velocities
    if(velocityDistributionType == "boltzmann") {
        generator.boltzmannDistributeVelocities(temperature, molecules);
    } else if(velocityDistributionType == "uniform") {
        generator.uniformDistributeVelocities(4, molecules);
    }
    // Set up force
    InteratomicForce* force = new InteratomicForce();
    force->setPotentialConstant(potentialConstant);
    force->setEnergyConstant(energyConstant);

    // Set up molecule system
    MoleculeSystem system;

    system.setInteratomicForce(force);

    // Set up file manager
    FileManager fileManager(&system);
    fileManager.setOutFileName(outFileName);
    fileManager.setUnitLength(unitLength);
    fileManager.setUnitTime(unitTime);
    fileManager.setUnitMass(unitMass);

    system.setFileManager(&fileManager);

    // Set up integrator
    Integrator* integrator = new VelocityVerletIntegrator(&system);
    integrator->setTimeStep(timeStep);
    system.setIntegrator(integrator);

    // Set up the rest of the system
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

