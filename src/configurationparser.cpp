#include "configurationparser.h"

// Local includes
#include <src/moleculesystem.h>
#include <src/generator.h>
#include <src/force/lennardjonesforce.h>
#include <src/force/vashishtatwoparticleforce.h>
#include <src/force/vashishtathreeparticleforce.h>
#include <src/integrator/integrator.h>
#include <src/integrator/velocityverletintegrator.h>
#include <src/integrator/eulercromerintegrator.h>
#include <src/filemanager.h>
#include <src/modifier/berendsenthermostat.h>
#include <src/modifier/andersenthermostat.h>
#include <src/force/constantforce.h>
#include <src/progressreporter.h>
#include <src/atom.h>

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
            try {
                config.readFile(configurationFileName.c_str());
            } catch (libconfig::FileIOException) {
                cerr << "Could not open config file " << configurationFileName << endl;
                exit(3121);
            }
        }
        world.barrier();
    }

    double unitLength = config.lookup("units.length");
    double unitTime = config.lookup("units.time");
    double unitMass = config.lookup("units.mass");
    double boltzmannConstant = config.lookup("units.boltzmannConstant");
    double unitVelocity = unitLength / unitTime;
    double unitTemperature = (unitVelocity * unitVelocity * unitMass) / (boltzmannConstant);
//    double unitTemperature = config.lookup("units.temperature");
    double unitForce = unitLength * unitMass / (unitTime * unitTime);
    double unitEnergy = unitLength * unitLength * unitMass / (unitTime * unitTime);

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
    energyConstant /= unitEnergy;

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

    unordered_map<int,AtomType> particleTypesByID;
    vector<AtomType> particleTypes;
    Setting& particleTypesSetting = config.lookup("particleTypes");
    int atomTypeIndex = 0;
    for(uint i = 0; true; i++) {
        try {
            AtomType atomType(atomTypeIndex);
            int id = -2;
            particleTypesSetting[i].lookupValue("id", id);
            atomType.setNumber(id);
            if(particleTypesSetting[i].exists("name")) {
                string name = "Noname";
                particleTypesSetting[i].lookupValue("name", name);
                atomType.setName(name);
            }
            if(particleTypesSetting[i].exists("abbreviation")) {
                string abbreviation = "NA";
                particleTypesSetting[i].lookupValue("abbreviation", abbreviation);
                atomType.setAbbreviation(abbreviation);
            }
            if(particleTypesSetting[i].exists("mass")) {
                double mass = 0;
                particleTypesSetting[i].lookupValue("mass", mass);
                mass /= unitMass;
                atomType.setMass(mass);
            }
            if(particleTypesSetting[i].exists("effectiveCharge")) {
                double effectiveCharge = 0;
                particleTypesSetting[i].lookupValue("effectiveCharge", effectiveCharge);
                atomType.setEffectiveCharge(effectiveCharge);
            }
            if(particleTypesSetting[i].exists("electronicPolarizability")) {
                double electronicPolarizability = 0;
                particleTypesSetting[i].lookupValue("electronicPolarizability", electronicPolarizability);
                atomType.setElectronicPolarizability(electronicPolarizability);
            }
            particleTypes.push_back(atomType);
            particleTypesByID[atomType.number()] = atomType;
            atomTypeIndex += 1;
        } catch(exception) {
            cout << "No more atom types found. Breaking." << endl;
            break;
        }
    }
    cout << "Particle types set:" << endl;
    for(AtomType type : particleTypes) {
        cout << type.name() << endl;
    }

    m_moleculeSystem->setParticleTypes(particleTypes);

    cout << "Loading initialization steps" << endl;
    Setting& initialization = config.lookup("initialization");
    for(uint i = 0; true; i++) {
        string initializationType;
        try {
            initialization[i].lookupValue("type", initializationType);
            cout << "Found initialization step " << initializationType << endl;
            if(initializationType == "fcc") {
                double b = initialization[i]["b"];
                b /= unitLength;
                int nCells = initialization[i]["nCells"];
                int particleType = initialization[i]["particleType"];
                AtomType atomType = particleTypesByID[particleType];
                vector<Atom*> atoms = generator.generateFcc(b, nCells, atomType);
                m_moleculeSystem->setBoundaries(generator.lastBoundaries());
                m_moleculeSystem->addAtoms(atoms);
                cout << "setbounds" << endl;
            } else if(initializationType == "moleculeTest") {
                m_moleculeSystem->setBoundaries(0,26,0,26,0,26);
                vector<Atom*> atoms;
                int idCounter = 1;
                int type = 0;
                for(int i = 0; i < 3; i++) {
                    for(int j = 0; j < 3; j++) {
                        type = ((i+j)%2);
                        for(int k = 0; k < 4; k++) {
                            Atom* atom;
                            if(type % 2) {
                                atom = new Atom(particleTypesByID[14]);
                            } else {
                                atom = new Atom(particleTypesByID[8]);
                            }
                            atom->setID(idCounter);
                            Vector3 position(i - 1,j,k);
                            position *= 3.0;
                            cout << position << endl;
                            atom->setPosition(position);
                            atoms.push_back(atom);
                            idCounter++;
                            type++;
                        }
                    }
                }
                m_moleculeSystem->addAtoms(atoms);
//                for(int i = 0; i < 3; i++) {
//                    for(int j = 0; j < 2; j++) {
//                        Atom* atom1;
//                        Atom* atom2;
//                        Atom* atom3;
//                        if(i % 2) {
//                            atom1 = new Atom(particleTypesByID[14]);
//                            atom2 = new Atom(particleTypesByID[14]);
//                            atom3 = new Atom(particleTypesByID[8 ]);
//                        } else {
//                            atom1 = new Atom(particleTypesByID[8]);
//                            atom2 = new Atom(particleTypesByID[8]);
//                            atom3 = new Atom(particleTypesByID[14]);
//                        }
//                        atom1->setID(1 + j * 3 + i * 2 * 3);
//                        atom2->setID(2 + j * 3 + i * 2 * 3);
//                        atom3->setID(3 + j * 3 + i * 2 * 3);
//                        double scale = 1.3;
//                        atom1->setPosition(Vector3(1.0, -0.2 + j*3, 1.0 + i * 2) * scale);
//                        atom2->setPosition(Vector3(3.0, -0.2 + j*3, 1.0 + i * 2) * scale);
//                        atom3->setPosition(Vector3(2.0, 1.5 + j*3, 1.0 + i * 2) * scale);
//    //                    atom1->setPosition(Vector3(1.0, 1.0, 1.0 + i * 2.5) * scale);
//    //                    atom2->setPosition(Vector3(3.0, 1.0, 1.0 + i * 2.5) * scale);
//    //                    atom3->setPosition(Vector3(2.0, 2.7, 1.0 + i * 2.5) * scale);
//                        vector<Atom*> atoms({atom1, atom2, atom3});
//                        m_moleculeSystem->addAtoms(atoms);
//                    }
//                }
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
        } catch(exception) {
            cout << "No more initialization steps found. Breaking." << endl;
            break;
        }
    }

    cout << "Setting up progress reporter" << endl;

    // Set up progressreporter
    ProgressReporter* reporter = new ProgressReporter("dragly", configurationFileName);
    m_moleculeSystem->setProgressReporter(reporter);

    // Set up force
//    LennardJonesForce* force = new LennardJonesForce();
//    force->setPotentialConstant(potentialConstant);
//    force->setEnergyConstant(energyConstant);
    VashishtaTwoParticleForce* twoParticleForce = new VashishtaTwoParticleForce();
    m_moleculeSystem->addTwoParticleForce(twoParticleForce);
    VashishtaThreeParticleForce* threeParticleForce = new VashishtaThreeParticleForce();
    m_moleculeSystem->addThreeParticleForce(threeParticleForce);

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
            BerendsenThermostat* thermostat = new BerendsenThermostat(m_moleculeSystem);
            double targetTemperature = modifiers[i]["targetTemperature"];
            targetTemperature /= unitTemperature;
            double relaxationTime = modifiers[i]["relaxationTime"];
            relaxationTime /= unitTime;
            thermostat->setTargetTemperature(targetTemperature);
            thermostat->setRelaxationTime(relaxationTime);
            m_moleculeSystem->addModifier(thermostat);
        } else if(modifierName == "andersenThermostat") {
            AndersenThermostat* thermostat = new AndersenThermostat(m_moleculeSystem);
            double targetTemperature = modifiers[i]["targetTemperature"];
            targetTemperature /= unitTemperature;
            double collisionTime = modifiers[i]["collisionTime"];
            collisionTime /= unitTime;
            thermostat->setTargetTemperature(targetTemperature);
            thermostat->setCollisionTime(collisionTime);
            m_moleculeSystem->addModifier(thermostat);
        } else if(modifierName == "constantForce") {
            ConstantForce* constantForce = new ConstantForce();
            Vector3 forceVector;
            for(int j = 0; j < 3; j++) {
                forceVector[j] = modifiers[i]["forceVector"][j];
            }
            forceVector /= unitForce;
            cout << "Setting constant force to " << forceVector << endl;
            constantForce->setForceVector(forceVector);
            m_moleculeSystem->addSingleParticleForce(constantForce);
        }
    }

    // Set up integrator
    Integrator* integrator = new VelocityVerletIntegrator(m_moleculeSystem);
    integrator->setTimeStep(timeStep);
    m_moleculeSystem->setIntegrator(integrator);

    // Set up the rest of the system
    //    m_moleculeSystem->loadConfiguration(&config); // TODO remove this
    cout << "addded" << endl;
//    m_moleculeSystem->setupCells(potentialConstant * 3);
    m_moleculeSystem->setupCells(4.43);
    m_moleculeSystem->setNSimulationSteps(nSimulationSteps);
    cout << "Setup cells" << endl;

    m_moleculeSystem->simulate();
}
