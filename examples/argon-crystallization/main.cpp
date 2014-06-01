#include "moleculesystem.h"
#include "moleculesystemcell.h"
#include "generator.h"
#include "atomtype.h"
#include "force/lennardjonesforce.h"
#include "force/fannthreeparticleforce.h"
#include "force/fanntwoparticleforce.h"
#include "modifier/andersenthermostat.h"
#include "modifier/berendsenthermostat.h"
#include "modifier/friction.h"
#include "atom.h"
#include "processor.h"
#include "utils/logging.h"
#include "force/threeparticleforce.h"
#include "integrator/velocityverletintegrator.h"
#include "filemanager.h"
#include "utils/setup.h"

#include <yaml-cpp/yaml.h>

using namespace std;

int main(int argc, char* argv[])
{    
    Setup setup(argc, argv);
    if(argc < 3) {
        cerr << "Too few arguments." << endl;
        cerr << "Usage: master-hydrogen-dissociation <config.yaml> <outfolder>" << endl;
        return 1;
    }

    ifstream fin(argv[1]);
    if(fin.fail()) {
        cerr << "Could not open the configuration file " << argv[1] << endl;
        return 1;
    }
    YAML::Parser parser(fin);

    double electronMassPerAtomicMass = 1822.8896;
    LOG(INFO) << "Example argon crystallization started";

    // Particle types
    AtomType argonType(0);
    argonType.setMass(39.948*electronMassPerAtomicMass);
    argonType.setNumber(18);
    argonType.setAbbreviation("Ar");

    vector<AtomType> particleTypes;
    particleTypes.push_back(argonType);

    double me = 9.10938291e-31; // kg
    double mp = 1.673823232368e-27; // kg
    double a0 = 5.2917721092e-11; // m
    double au_rho = me / (a0*a0*a0);
    double au_T = 3.1577464e5; // K
    double au_t = 2.418884326505e-17; // s
    double au_E = 4.35974417e-18; // J

    double sigma = 3.4e-10 / a0;
    double epsilon = 1.65e-21 / au_E;

    YAML::Node rootNode;
    parser.GetNextDocument(rootNode);
    double cellCutoff = 2.3*sigma;
    double cutoffRadius = 2.3*sigma;
    double targetTemperature = 0.0;
    double boilingTemperature = 0.0;
    int timeStepsInRun = 100000;
    double density = 0.0;
    rootNode["temperature"] >> targetTemperature;
    targetTemperature /= au_T;
    rootNode["boiling_temperature"] >> boilingTemperature;
    boilingTemperature /= au_T;
    double timeStep = 10.0;

    // Read optional params or catch exception
    try {
        rootNode["cell_cutoff"] >> cellCutoff;
    } catch(YAML::TypedKeyNotFound<std::string>& ) {
    }
    try {
        rootNode["cell_cutoff"] >> cellCutoff;
    } catch(YAML::TypedKeyNotFound<std::string>& ) {
    }
    try {
        rootNode["timestep_count"] >> timeStepsInRun;
    } catch(YAML::TypedKeyNotFound<std::string>& ) {
    }
    try {
        rootNode["cutoff_radius"] >> cutoffRadius;
    } catch(YAML::TypedKeyNotFound<std::string>& ) {
    }
    try {
        rootNode["timestep"] >> timeStep;
    } catch(YAML::TypedKeyNotFound<std::string>& ) {
    }

    cout << "Temperature: " << targetTemperature << endl;
    cout << "Density: " << density << endl;
    cout << "cellCutoff: " << cellCutoff << endl;
    cout << "cutoffRadius: " << cutoffRadius << endl;
    cout << "timeStep: " << timeStep << endl;

//    double sideLength = 31.96e-10 / a0;

    int nX = 16;
    int nY = nX;
    int nZ = nX;
    double sideLength = 31.96e-10 / a0 * nX / 8.0;

    cout << "nX: " << nX << endl;
    cout << "nY: " << nY << endl;
    cout << "nZ: " << nZ << endl;

    double spacingx = sideLength / nX;
    double spacingy = sideLength / nY;
    double spacingz = sideLength / nZ;

    Generator generator;
    int idCounter = 0;
    vector<Atom*> atoms;
    for(int i = 0; i < nZ; i++) {
        for (int j = 0; j < nY; j++) {
            for (int k = 0; k < nX; ++k) {
                Atom *argonAtom1 = new Atom(argonType);
                argonAtom1->setPosition(Vector3(0.5 + spacingx*k, 0.5 + spacingy*j, 0.5 + spacingz*i));
                argonAtom1->setID(1 + idCounter);
                atoms.push_back(argonAtom1);
                idCounter += 1;
            }
        }
    }
    cout << "Number of atoms: " << atoms.size() << endl;

    generator.boltzmannDistributeVelocities(1e-9, atoms);

    MoleculeSystem system;
    system.setPeriodicity(true, true, true);
    system.setBoundaries(0, sideLength);
    system.setParticleTypes(particleTypes);
    system.addAtoms(atoms);

    LennardJonesForce lennardJonesForce;
    lennardJonesForce.setCalculatePotentialEnabled(true);
    lennardJonesForce.setCalculatePressureEnabled(true);
    lennardJonesForce.setCutoffRadius(2.3*sigma);
    lennardJonesForce.setPotentialConstant(sigma);
    lennardJonesForce.setEnergyConstant(epsilon);
    system.setTwoParticleForce(&lennardJonesForce);

    VelocityVerletIntegrator integrator(&system);
    integrator.setTimeStep(timeStep);
    system.setIntegrator(&integrator);

    BerendsenThermostat thermostat(&system);
    system.addModifier(&thermostat);

    FileManager fileManager(&system);
    fileManager.setOutFileName(argv[2]);
    fileManager.setUnitLength(a0);
    fileManager.setUnitMass(me);
    fileManager.setUnitTemperature(au_T);
    fileManager.setUnitTime(au_t);
    system.setFileManager(&fileManager);

    system.setSaveEnabled(true);
    system.setOutputEnabled(true);

    system.setupCells(cellCutoff);
    system.save();

    cout << "Cells: " << system.nCells() << endl;

    system.setSaveEveryNSteps(100);

    // Boiling
    cout << "Boiling to " << boilingTemperature << endl;
    thermostat.setTargetTemperature(boilingTemperature);
    thermostat.setRelaxationTime(100 * timeStep);
    system.setNSimulationSteps(10000);
    system.simulate();
    generator.removeLinearMomentum(system.atoms());

    // Final thermalization
    cout << "Final thermalization to " << targetTemperature << endl;
    int thermalStepCount = 1000;
    for(int i = 0; i < thermalStepCount; i++) {
        double stepRatio = double(i) / double(thermalStepCount);
        thermostat.setTargetTemperature(boilingTemperature - (boilingTemperature - targetTemperature) * stepRatio);
        thermostat.setRelaxationTime(1e-13 / au_t);
        cout << "Thermalizing to " << thermostat.targetTemperature() << endl;
        system.setNSimulationSteps(timeStepsInRun / thermalStepCount);
        system.simulate();
        generator.removeLinearMomentum(system.atoms());

        system.removeModifier(&thermostat);
        system.setNSimulationSteps(100);
        system.simulate();
        system.addModifier(&thermostat);
    }

    LOG(INFO) << "FannForceWater test complete";
    return 0;
}

