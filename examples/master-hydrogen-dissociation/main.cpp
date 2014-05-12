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
    LOG(INFO) << "Example master hydrogen dissociation started";

    // Particle types
    AtomType hydrogenType(0);
    hydrogenType.setMass(1.008*electronMassPerAtomicMass);
    hydrogenType.setNumber(1);
    hydrogenType.setAbbreviation("H");

    vector<AtomType> particleTypes;
    particleTypes.push_back(hydrogenType);

    vector<Atom *> atoms;

    double me = 9.10938291e-31; // kg
    double mp = 1.673823232368e-27; // kg
    double a0 = 5.2917721092e-11; // m
    double au_rho = me / (a0*a0*a0);
    double au_T = 3.1577464e5; // K
    double au_t = 2.418884326505e-17; // s

    YAML::Node rootNode;
    parser.GetNextDocument(rootNode);
    double cellCutoff = 12.0;
    double cutoffRadius = 12.0;
    double targetTemperature = 0.0;
    double density = 0.0;
    rootNode["temperature"] >> targetTemperature;
    targetTemperature /= au_T;
    rootNode["density"] >> density;
    density /= au_rho;
    double timeStep = 10.0;

    // Read optional params or catch exception
    try {
        rootNode["cell_cutoff"] >> cellCutoff;
    } catch(YAML::TypedKeyNotFound<std::string>& ) {
    }
    try {
        rootNode["cutoff_radius"] >> cutoffRadius;
    } catch(YAML::TypedKeyNotFound<std::string>& ) {
    }
    try {
        rootNode["time_step"] >> timeStep;
    } catch(YAML::TypedKeyNotFound<std::string>& ) {
    }

    double numberDensity = density / (mp/me);
    cout << "Temperature: " << targetTemperature << endl;
    cout << "Density: " << density << endl;
    cout << "cellCutoff: " << cellCutoff << endl;
    cout << "cutoffRadius: " << cutoffRadius << endl;
    cout << "timeStep: " << timeStep << endl;

    double initialTemperatureBoiling = 15600 / au_T;
    double initialTemperatureLiquid = 300 / au_T;

    int atomCountIn = 0;
    rootNode["atom_count"] >> atomCountIn;

    double volume = atomCountIn / numberDensity;
    double sideLengthX = pow(4*volume, 1./3.);
    double sideLengthY = sideLengthX / 2.0;
    double sideLengthZ = sideLengthX / 2.0;

    int nX = round(pow(4*atomCountIn / 2.0, 1./3.));
    int nY = round(nX / 2.0);
    int nZ = round(nX / 2.0);

    cout << "nX: " << nX << endl;
    cout << "nY: " << nY << endl;
    cout << "nZ: " << nZ << endl;

    double spacingx = sideLengthX / nX;
    double spacingy = sideLengthY / nY;
    double spacingz = sideLengthZ / nZ;
    cout << "Side lengths: " << sideLengthX << " " << sideLengthY << " " << sideLengthZ << endl;

    int idCounter = 0;
    for(int i = 0; i < nZ; i++) {
        for (int j = 0; j < nY; j++) {
            for (int k = 0; k < nX; ++k) {
                if(idCounter >= atomCountIn) {
                    continue;
                }
                Atom *hydrogenAtom1 = new Atom(hydrogenType);
                hydrogenAtom1->setPosition(Vector3(0.5 + spacingx*k, 0.5 + spacingy*j, 0.5 + spacingz*i));
                hydrogenAtom1->setID(1 + idCounter);
                atoms.push_back(hydrogenAtom1);
                Atom *hydrogenAtom2 = new Atom(hydrogenType);
                hydrogenAtom2->setPosition(Vector3(0.5 + spacingx*k + 1.4, 0.5 + spacingy*j, 0.5 + spacingz*i));
                hydrogenAtom2->setID(2 + idCounter);
                atoms.push_back(hydrogenAtom2);
                idCounter += 2;
            }
        }
    }
    cout << "Number of atoms: " << atoms.size() << endl;
    cout << "Volume: " << volume << endl;
    double totalMass = 0.0;
    for(Atom* atom : atoms) {
        totalMass += atom->mass();
    }
    cout << "Total mass: " << totalMass << endl;
    cout << "Number density check: " << atoms.size() / (sideLengthX*sideLengthY*sideLengthZ) << endl;
    for(Atom* atom : atoms) {
        atom->setPosition(atom->position() + Vector3(sideLengthX / 2, sideLengthY / 2, sideLengthZ / 2));
    }

    Generator gen;
    gen.boltzmannDistributeVelocities(1e-9, atoms);

    MoleculeSystem system;
    system.setPeriodicity(true, true, true);
    system.setBoundaries(0.0, sideLengthX, 0.0, sideLengthY, 0.0, sideLengthZ);
    system.setParticleTypes(particleTypes);
    system.addAtoms(atoms);

    FannTwoParticleForce testForce2;
    testForce2.setCutoffRadius(cutoffRadius);
    testForce2.addNetwork(hydrogenType, hydrogenType,
                          "/home/svenni/Dropbox/studies/master/results/fann_train/20140507-204601/fann_network.net",
                          "/home/svenni/Dropbox/studies/master/results/fann_train/20140507-204601/bounds.fann");
    testForce2.addNetwork(hydrogenType, hydrogenType,
                          "/home/svenni/Dropbox/studies/master/results/fann_train/20140507-200839/fann_network.net",
                          "/home/svenni/Dropbox/studies/master/results/fann_train/20140507-200839/bounds.fann");
    system.setTwoParticleForce(&testForce2);

    FannThreeParticleForce testForce3;
    testForce3.setCutoffRadius(fmin(cellCutoff / 2.0, cutoffRadius));
    testForce3.loadNetwork("/home/svenni/Dropbox/studies/master/results/fann_train/20140512-200948/fann_network.net",
                           "/home/svenni/Dropbox/studies/master/results/fann_train/20140512-200948/bounds.fann",
                           5.0);
    testForce3.loadNetwork("/home/svenni/Dropbox/studies/master/results/fann_train/20140512-182447/fann_network.net",
                           "/home/svenni/Dropbox/studies/master/results/fann_train/20140512-182447/bounds.fann");
    system.setThreeParticleForce(&testForce3);

    VelocityVerletIntegrator integrator(&system);
    integrator.setTimeStep(timeStep);
    system.setIntegrator(&integrator);

    BerendsenThermostat thermostat(&system);

    FileManager fileManager(&system);
    fileManager.setOutFileName(argv[2]);
    fileManager.setUnitLength(a0);
    fileManager.setUnitMass(mp);
    fileManager.setUnitTemperature(au_T);
    fileManager.setUnitTime(au_t);
    system.setFileManager(&fileManager);

    system.setSaveEnabled(true);
    system.setOutputEnabled(true);

    system.setupCells(cellCutoff);
    system.save();

    cout << "Cells: " << system.nCells() << endl;

    int timeStepsInRun = 1000000;

    system.setSaveEveryNSteps(10000);
    // Equilibrate initial
    system.addModifier(&thermostat);
    thermostat.setTargetTemperature(1e-3);
    thermostat.setRelaxationTime(timeStep * 10);
    system.setNSimulationSteps(timeStepsInRun / 50);
    system.simulate();

    // Cooking
    integrator.setTimeStep(timeStep*0.1);
    thermostat.setTargetTemperature(initialTemperatureBoiling);
    thermostat.setRelaxationTime(timeStep * 10.0);
    system.setNSimulationSteps(timeStepsInRun / 50);
    system.simulate();
    gen.removeLinearMomentum(system.atoms());
    integrator.setTimeStep(timeStep);

    if(initialTemperatureLiquid > targetTemperature) {
        // Liquidizing
        thermostat.setTargetTemperature(initialTemperatureLiquid);
        thermostat.setRelaxationTime(timeStepsInRun / 20);
        system.setNSimulationSteps(timeStepsInRun / 10);
        system.simulate();
        gen.removeLinearMomentum(system.atoms());
    }

    // Equilibration
    system.removeModifier(&thermostat);
    system.setNSimulationSteps(timeStepsInRun / 10);
    system.simulate();

    // Final thermalization
    system.addModifier(&thermostat);
    thermostat.setTargetTemperature(targetTemperature);
    thermostat.setRelaxationTime(timeStepsInRun / 20);
    system.setNSimulationSteps(timeStepsInRun / 10);
    system.simulate();
    gen.removeLinearMomentum(system.atoms());

    // Final run
    system.setSaveEveryNSteps(1000);
    system.setNSimulationSteps(timeStepsInRun);
    system.simulate();

    // Final dynamics
    system.setSaveEveryNSteps(1);
    system.setNSimulationSteps(timeStepsInRun / 1000);
    system.simulate();

    //        CHECK_EQUAL(3, system.nAtomsTotal());

    LOG(INFO) << "FannForceWater test complete";
    return 0;
}

