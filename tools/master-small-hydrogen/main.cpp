#include <unittest++/UnitTest++.h>

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

#ifdef MD_USE_MPI
#include <boost/mpi.hpp>
namespace mpi = boost::mpi;
#endif
double electronMassPerAtomicMass = 1822.8896;

int main(int argc, char* argv[]) {
    setlocale(LC_ALL, "C");
#ifdef MD_USE_GLOG
    FLAGS_log_dir = ".";
    google::InitGoogleLogging(argv[0]);
#endif
#ifdef MD_USE_MPI
    mpi::environment env(argc, argv);
    mpi::communicator world;
    (void)env;
    (void)world;
#endif

    LOG(INFO) << "Master small hydrogen started";

    // Particle types
    AtomType hydrogenType(0);
    hydrogenType.setMass(1.008*electronMassPerAtomicMass);
    hydrogenType.setNumber(1);
    hydrogenType.setAbbreviation("H");

    vector<AtomType> particleTypes;
    particleTypes.push_back(hydrogenType);

    vector<Atom *> atoms;

    bool friction = false;
    bool startVelocities = true;
    bool thermo = true;
    bool periodic = true;
    double cellCutoff = 12.0;
    double cutoffRadius = 12.0;
    double sideLength = 20.0;
    double targetTemperature = 10.0e-5;

    if(!periodic) {
        sideLength = 40;
    }
    double numberDensity = 0.00627; // 70.82 kg / m^3
    double numberDensityHalf = numberDensity / 2.0;
    //                targetTemperature = 3.16681542254e-05; // 10 K (solid)
    //            targetTemperature = 4.43354159156e-05; // 14 K (phase change)
    //    targetTemperature = 5.0e-5; // 17 K (liquid)
    //            targetTemperature = 0.000158340771127; // 50 K (liquid)
    targetTemperature = 0.000475022313381; // 150 K (liquid)
    //            targetTemperature = 0.000494;

    int nx = 3;
    int ny = 2;
    int nz = 2;

    double volume = (nx * ny * nz) / numberDensityHalf;

    sideLength =  pow(volume, 1.0/3.0);

    cout << "Side length: " << sideLength << endl;
    double spacingx = sideLength / nx;
    double spacingy = sideLength / ny;
    double spacingz = sideLength / nz;
    if(!periodic) {
        spacingx = sideLength / 5 / nx;
        spacingy = sideLength / 5 / ny;
        spacingz = sideLength / 5 / nz;
    }

    int idCounter = 0;
    for(int i = 0; i < nz; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nx; ++k) {
                Atom *hydrogenAtom1 = new Atom(hydrogenType);
                hydrogenAtom1->setPosition(Vector3(0.5 + spacingx*k, 0.5 + spacingy*j, 0.5 + spacingz*i));
                hydrogenAtom1->setID(1 + idCounter);
                atoms.push_back(hydrogenAtom1);
                idCounter += 1;
            }
        }
    }
    cout << "Number of atoms: " << atoms.size() << endl;
    cout << "Volume: " << sideLength*sideLength*sideLength << endl;
    double totalMass = 0.0;
    for(Atom* atom : atoms) {
        totalMass += atom->mass();
    }
    cout << "Total mass: " << totalMass << endl;
    cout << "Density: " << totalMass / (sideLength*sideLength*sideLength) << endl;
    if(!periodic) {
        for(Atom* atom : atoms) {
            atom->setPosition(atom->position() + Vector3(sideLength / 2, sideLength / 2, sideLength / 2));
        }
    }


    if(startVelocities) {
        Generator gen;
        if(friction) {
            gen.boltzmannDistributeVelocities(targetTemperature*10, atoms);
        } else {
            gen.boltzmannDistributeVelocities(targetTemperature, atoms);
        }
    }

    MoleculeSystem system;

    FileManager fileManager(&system);
    fileManager.setOutFileName("/tmp/fannforce/hydrogen101/atoms*.bin");
    fileManager.setUnitLength(5.2917721092e-11);
    fileManager.setUnitMass(9.10938291e-31);
    fileManager.setUnitTime(2.418884326505e-17);
    fileManager.setUnitTemperature(3.1577464e5);

    system.setFileManager(&fileManager);
    system.load("/home/svenni/Dropbox/compphys/people/Milad/hydrogen-comparison/hydrogen-liquid-12-atoms/atoms140000.bin");
    system.setPeriodicity(periodic, periodic, periodic);
    system.setBoundaries(0.0, sideLength, 0.0, sideLength, 0.0, sideLength);
    system.setParticleTypes(particleTypes);

    // TODO Load atom types properly
    for(Atom* atom : system.atoms()) {
        atom->setAtomType(hydrogenType);
        atom->clearForcePotentialPressure();
    }

    FannTwoParticleForce testForce2;
    testForce2.setCutoffRadius(cutoffRadius);
    //        testForce2.addNetwork(hydrogenType, hydrogenType,
    //                              "/home/svenni/Dropbox/studies/master/results/fann_train/20140427-141531/fann_network.net",
    //                              "/home/svenni/Dropbox/studies/master/results/fann_train/20140427-141531/bounds.fann");
    testForce2.addNetwork(hydrogenType, hydrogenType,
                          "/home/svenni/Dropbox/studies/master/results/fann_train/20140507-200839/fann_network.net",
                          "/home/svenni/Dropbox/studies/master/results/fann_train/20140507-200839/bounds.fann");
    testForce2.addNetwork(hydrogenType, hydrogenType,
                          "/home/svenni/Dropbox/studies/master/results/fann_train/20140507-204601/fann_network.net",
                          "/home/svenni/Dropbox/studies/master/results/fann_train/20140507-204601/bounds.fann");

    //        testForce2.addNetwork(hydrogenType, hydrogenType,
    //                              "/home/svenni/Dropbox/projects/programming/fann-md/fann-md/tools/train/tmp/two/fann_network.net",
    //                              "/home/svenni/Dropbox/projects/programming/fann-md/fann-md/tools/train/tmp/two/bounds.fann");


    system.setTwoParticleForce(&testForce2);

    FannThreeParticleForce testForce3;
    testForce3.setCutoffRadius(fmin(cellCutoff / 2.0, cutoffRadius));
    testForce3.loadNetwork("/home/svenni/Dropbox/studies/master/results/fann_train/20140428-194643/fann_network_0.net",
                           "/home/svenni/Dropbox/studies/master/results/fann_train/20140428-194643/bounds.fann");

    system.setThreeParticleForce(&testForce3);

    VelocityVerletIntegrator integrator(&system);
    // Minimum for energy conservation: 0.001
    integrator.setTimeStep(10.0);
    system.setIntegrator(&integrator);

    BerendsenThermostat thermostat(&system);
    if(thermo) {
        thermostat.setRelaxationTime(1000.0);
    }

    Friction frictionModifier(&system);

    system.setSaveEnabled(true);
    system.setSaveEveryNSteps(1000);
    system.setOutputEnabled(true);

    system.setupCells(cellCutoff);

    system.save();


    cout << "Cells: " << system.nCells();

//    system.addModifier(&thermostat);
//    thermostat.setTargetTemperature(100*targetTemperature);
//    system.setNSimulationSteps(40000);
//    system.simulate();

//    system.removeModifier(&thermostat);
//    system.setNSimulationSteps(10000);
//    system.simulate();

//    system.addModifier(&thermostat);
//    thermostat.setTargetTemperature(targetTemperature);
    system.setNSimulationSteps(1000000);
    system.simulate();

//    system.setNSimulationSteps(100000);
//    system.simulate();

    LOG(INFO) << "FannForceWater test complete";
}
