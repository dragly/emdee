#include <src/moleculesystem.h>
#include <src/generator.h>
#include <src/interatomicforce.h>
#include <src/integrator/integrator.h>
#include <src/integrator/velocityverletintegrator.h>
#include <src/integrator/eulercromerintegrator.h>
#include <src/filemanager.h>
#include <src/modifier/berendsenthermostat.h>
#include <src/modifier/andersenthermostat.h>

#include <iostream>
#include <libconfig.h++>

using namespace std;
using namespace libconfig;

int main(int argc, char** argv)
{
    Config config;
    string configFileName = "testconfig.cfg";
    if(argc > 1) {
        configFileName = argv[1];
    }
    config.readFile(configFileName.c_str());

    double unitLength = config.lookup("units.length");
    double unitTime = config.lookup("units.time");
    double unitMass = config.lookup("units.mass");
    double unitTemperature = config.lookup("units.temperature");
//    double boltzmannConstant = config.lookup("units.boltzmannConstant");
//    boltzmannConstant /= (unitLength*unitLength * unitMass / (unitTime*unitTime));

    // TODO Why must the boltzmannConstant be set to one here?
//    boltzmannConstant = 1;

    // File manager config
    string outFileName;
    config.lookupValue("simulation.saveFileName", outFileName);

    int nSimulationSteps = config.lookup("simulation.nSimulationSteps");
    double timeStep = config.lookup("integrator.timeStep");
    timeStep /= unitTime;

    Generator generator;
//    generator.setBoltzmannConstant(boltzmannConstant);
//    generator.loadConfiguration(&config);

    // Generator specific config
//    string velocityDistributionType;
//    config.lookupValue("generator.velocity.distribution", velocityDistributionType);
//    int nCells = config.lookup("generator.fcc.nCells");
//    double b = config.lookup("generator.fcc.b");
//    double temperature = config.lookup("system.initialTemperature");
//    temperature /= unitTemperature;
//    b /= unitLength;

    // Force and potentials
    double potentialConstant = config.lookup("system.potentialConstant");
    potentialConstant /= unitLength;
    double energyConstant = config.lookup("system.energyConstant");
    energyConstant /= (unitMass * unitLength * unitLength / (unitTime * unitTime));

    // Set up initial system
    // Set up molecule system
    MoleculeSystem system;

    // Set up file manager
    FileManager fileManager(&system);
    fileManager.setOutFileName(outFileName);
    fileManager.setUnitLength(unitLength);
    fileManager.setUnitTime(unitTime);
    fileManager.setUnitMass(unitMass);
    fileManager.setUnitTemperature(unitTemperature);

    system.setFileManager(&fileManager);

    Setting& initialization = config.lookup("initialization");
    for(uint i = 0; true; i++) {
        string initializationType;
        try {
            initialization[i].lookupValue("type", initializationType);
        } catch(exception) {
            cout << "No more initialization steps found. Breaking." << endl;
            break;
        }
        cout << "Found initialization step " << initializationType << endl;
        if(initializationType == "fcc") {
            double b = initialization[i]["b"];
            b /= unitLength;
            int nCells = initialization[i]["nCells"];
            vector<Molecule*> molecules = generator.generateFcc(b, nCells, AtomType::argon());
            system.addMolecules(molecules);
            system.setBoundaries(generator.lastBoundaries());
            cout << "setbounds" << endl;
        } else if(initializationType == "boltzmannVelocity") {
            double initialTemperature = initialization[i]["initialTemperature"];
            initialTemperature /= unitTemperature;
            generator.boltzmannDistributeVelocities(initialTemperature, system.molecules());
        } else if(initializationType == "uniformVelocity") {
            double maxVelocity = initialization[i]["maxVelocity"];
            maxVelocity /= (unitLength / unitTime);
            generator.uniformDistributeVelocities(maxVelocity, system.molecules());
        } else if(initializationType == "loadFile") {
            string fileName = initialization[i]["fileName"];
            if(!system.load(fileName)) {
                cerr << "Could not load file " << fileName << endl;
                throw(new exception);
            }
        }
    }
//    vector<Molecule*> molecules = generator.generateFcc(b, nCells, AtomType::argon());

    // Set up force
    InteratomicForce* force = new InteratomicForce();
    force->setPotentialConstant(potentialConstant);
    force->setEnergyConstant(energyConstant);

    system.setInteratomicForce(force);
//    system.setBoltzmannConstant(boltzmannConstant);

    // Set up modifiers
    Setting& modifiers = config.lookup("modifiers");
    for(uint i = 0; true; i++) {
        string modifierName;
        try {
            modifiers[i].lookupValue("type", modifierName);
        } catch(exception) {
            cout << "No more modifiers found. Breaking." << endl;
            break;
        }
        cout << "Found modifier " << modifierName << endl;
        if(modifierName == "berendsenThermostat") {
            BerendsenThermostat thermostat(&system);
            double targetTemperature = modifiers[i]["targetTemperature"];
            targetTemperature /= unitTemperature;
            double relaxationTime = modifiers[i]["relaxationTime"];
            relaxationTime /= unitTime;
            thermostat.setTargetTemperature(targetTemperature);
            thermostat.setRelaxationTime(relaxationTime);
            system.addModifier(&thermostat);
        } else if(modifierName == "andersenThermostat") {
            AndersenThermostat thermostat(&system);
            double targetTemperature = modifiers[i]["targetTemperature"];
            targetTemperature /= unitTemperature;
            double collisionTime = modifiers[i]["collisionTime"];
            collisionTime /= unitTime;
            thermostat.setTargetTemperature(targetTemperature);
            thermostat.setCollisionTime(collisionTime);
            system.addModifier(&thermostat);
        }
    }

    // Set up integrator
    Integrator* integrator = new VelocityVerletIntegrator(&system);
    integrator->setTimeStep(timeStep);
    system.setIntegrator(integrator);

    // Set up the rest of the system
    system.loadConfiguration(&config); // TODO remove this
    cout << "addded" << endl;
    system.setupCells(potentialConstant * 3);
    cout << "Setup cells" << endl;
    system.simulate(nSimulationSteps);
    ofstream test;
    test.open("test");
    return 0;
}

