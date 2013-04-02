#include "configurationparser.h"

// Local includes
#include <src/moleculesystem.h>
#include <src/generator.h>
#include <src/force/lennardjonesforce.h>
#include <src/integrator/integrator.h>
#include <src/integrator/velocityverletintegrator.h>
#include <src/integrator/eulercromerintegrator.h>
#include <src/filemanager.h>
#include <src/modifier/berendsenthermostat.h>
#include <src/modifier/andersenthermostat.h>

// System includes
#include <boost/mpi.hpp>
namespace mpi = boost::mpi;

#include <string>
using namespace std;

ConfigurationParser::ConfigurationParser(MoleculeSystem *moleculeSystem) :
    m_moleculeSystem(moleculeSystem)
{
}

void ConfigurationParser::runConfiguration(string configurationFileName) {
    Config config;
    // Make sure only one processor reads the config file at a time
    mpi::communicator world;
    for(int rank = 0; rank < world.size(); rank++) {
        if(world.rank() == rank) {
            config.readFile(configurationFileName.c_str());
        }
        world.barrier();
    }

    double unitLength = config.lookup("units.length");
    double unitTime = config.lookup("units.time");
    double unitMass = config.lookup("units.mass");
    double unitTemperature = config.lookup("units.temperature");

    // File manager config
    string outFileName;
    config.lookupValue("simulation.saveFileName", outFileName);
    bool isSaveEnabled = config.lookup("simulation.saveEnabled");
    bool isCalcualatePressureEnabled = config.lookup("simulation.calculatePressure");
    bool isCalcualatePotentialEnabled = config.lookup("simulation.calculatePotential");
    int saveEveryNSteps = config.lookup("simulation.saveEveryNSteps");

    int nSimulationSteps = config.lookup("simulation.nSimulationSteps");
    double timeStep = config.lookup("integrator.timeStep");
    timeStep /= unitTime;

    Generator generator;

    // Force and potentials
    double potentialConstant = config.lookup("system.potentialConstant");
    potentialConstant /= unitLength;
    double energyConstant = config.lookup("system.energyConstant");
    energyConstant /= (unitMass * unitLength * unitLength / (unitTime * unitTime));

    // Set up initial system
    // Set up molecule system

    // Set up file manager
    FileManager fileManager(m_moleculeSystem);
    fileManager.setOutFileName(outFileName);
    fileManager.setUnitLength(unitLength);
    fileManager.setUnitTime(unitTime);
    fileManager.setUnitMass(unitMass);
    fileManager.setUnitTemperature(unitTemperature);
    string configurationName = configurationFileName.substr(configurationFileName.find_last_of("/") + 1);
    fileManager.setConfigurationName(configurationName);

    m_moleculeSystem->setFileManager(&fileManager);
    m_moleculeSystem->setSaveEnabled(isSaveEnabled);
    m_moleculeSystem->setCalculatePotentialEnabled(isCalcualatePotentialEnabled);
    m_moleculeSystem->setCalculatePressureEnabled(isCalcualatePressureEnabled);
    m_moleculeSystem->setSaveEveryNSteps(saveEveryNSteps);

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
            vector<Atom*> atoms = generator.generateFcc(b, nCells, AtomType::argon());
            m_moleculeSystem->setBoundaries(generator.lastBoundaries());
            m_moleculeSystem->addAtoms(atoms);
            cout << "setbounds" << endl;
        } else if(initializationType == "boltzmannVelocity") {
            double initialTemperature = initialization[i]["initialTemperature"];
            initialTemperature /= unitTemperature;
            generator.boltzmannDistributeVelocities(initialTemperature, m_moleculeSystem->atoms());
        } else if(initializationType == "uniformVelocity") {
            double maxVelocity = initialization[i]["maxVelocity"];
            maxVelocity /= (unitLength / unitTime);
            generator.uniformDistributeVelocities(maxVelocity, m_moleculeSystem->atoms());
        } else if(initializationType == "loadFile") {
            string fileName = initialization[i]["fileName"];
            if(!m_moleculeSystem->load(fileName)) {
                cerr << "Could not load file " << fileName << endl;
                throw(new exception);
            }
        }
    }

    // Set up force
    LennardJonesForce* force = new LennardJonesForce();
    force->setPotentialConstant(potentialConstant);
    force->setEnergyConstant(energyConstant);

    m_moleculeSystem->setInteratomicForce(force);

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
            BerendsenThermostat thermostat(m_moleculeSystem);
            double targetTemperature = modifiers[i]["targetTemperature"];
            targetTemperature /= unitTemperature;
            double relaxationTime = modifiers[i]["relaxationTime"];
            relaxationTime /= unitTime;
            thermostat.setTargetTemperature(targetTemperature);
            thermostat.setRelaxationTime(relaxationTime);
            m_moleculeSystem->addModifier(&thermostat);
        } else if(modifierName == "andersenThermostat") {
            AndersenThermostat thermostat(m_moleculeSystem);
            double targetTemperature = modifiers[i]["targetTemperature"];
            targetTemperature /= unitTemperature;
            double collisionTime = modifiers[i]["collisionTime"];
            collisionTime /= unitTime;
            thermostat.setTargetTemperature(targetTemperature);
            thermostat.setCollisionTime(collisionTime);
            m_moleculeSystem->addModifier(&thermostat);
        }
    }

    // Set up integrator
    Integrator* integrator = new VelocityVerletIntegrator(m_moleculeSystem);
    integrator->setTimeStep(timeStep);
    m_moleculeSystem->setIntegrator(integrator);

    // Set up the rest of the system
//    m_moleculeSystem->loadConfiguration(&config); // TODO remove this
    cout << "addded" << endl;
    m_moleculeSystem->setupCells(potentialConstant * 3);
    m_moleculeSystem->setNSimulationSteps(nSimulationSteps);
    cout << "Setup cells" << endl;

    m_moleculeSystem->simulate();
}
