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

double electronMassPerAtomicMass = 1822.8896;

SUITE(FannForceSystem) {
    TEST(Dummy) {

    }
    TEST(FannForceHydrogen)
    {
        LOG(INFO) << "FannForceHydrogen test started";

        // Particle types
        AtomType hydrogenType(0);
        hydrogenType.setMass(1.008*electronMassPerAtomicMass);
        hydrogenType.setNumber(1);
        hydrogenType.setAbbreviation("H");

        vector<AtomType> particleTypes;
        particleTypes.push_back(hydrogenType);

        vector<Atom *> atoms;

        bool friction = false;
        bool startVelocities = false;
        bool thermo = true;
        bool periodic = false;
        double cellCutoff = 100.0;
        double cutoffRadius = 12.0;
        double sideLength = 20.0;
        double targetTemperature = 4e-5;

        if(!periodic) {
            sideLength = 40;
        }

        periodic = true;
        friction = false;
        thermo = true;
        startVelocities = true;
        //            double numberDensity = 0.00627; // 70.82 kg / m^3
        double numberDensity = 0.00169; // rho* = 0.0011
        //            double numberDensity = 0.00175;
        //            double numberDensity = 7.95918040174197e-06; // 0.08988 kg / m3
        double numberDensityHalf = numberDensity / 2.0;
        //            targetTemperature = 5e-06; // (solid)
        //            targetTemperature = 3.16681542254e-05; // 10 K (solid)
//        targetTemperature = 4.43354159156e-05; // 14 K (phase change)
        //            targetTemperature = 5.0e-5; // 17 K (liquid)
        //            targetTemperature = 0.000158340771127; // 50 K (liquid)
//        targetTemperature = 0.000475022313381; // 150 K (liquid)
//        targetTemperature = 0.0475022313381; // 15000 K (dense liquid?)
        //            targetTemperature = 0.000494;

        int nx = 3;
        int ny = nx;
        int nz = nx;

        double volume = 10000;

        cout << "Side length: " << sideLength << endl;
        double spacingx = sideLength / nx;
        double spacingy = sideLength / ny;
        double spacingz = sideLength / nz;
        if(!periodic) {
            spacingx = 10;
            spacingy = 10;
            spacingz = 10;
        }

        int idCounter = 0;
        for(int i = 0; i < nz; i++) {
            for (int j = 0; j < ny; j++) {
                for (int k = 0; k < nx; ++k) {
                    Atom *hydrogenAtom1 = new Atom(hydrogenType);
                    hydrogenAtom1->setPosition(Vector3(0.5 + spacingx*k, 0.5 + spacingy*j, 0.5 + spacingz*i));
                    hydrogenAtom1->setID(1 + idCounter);
                    atoms.push_back(hydrogenAtom1);
                    Atom *hydrogenAtom2 = new Atom(hydrogenType);
                    hydrogenAtom2->setPosition(Vector3(0.5 + spacingx*k + 1.4, 0.5 + spacingy*j, 0.5 + spacingz*i));
                    hydrogenAtom2->setID(2 + idCounter);
                    atoms.push_back(hydrogenAtom2);
                    idCounter += 2;
                    cout << hydrogenAtom1->position() << endl;
                }
            }
        }
        for(Atom* atom : atoms) {
            atom->setPosition(atom->position() + Vector3(sideLength / 2, sideLength / 2, sideLength / 2));
        }

        Generator gen;
        gen.boltzmannDistributeVelocities(1e-6, atoms);

        MoleculeSystem system;
        system.setPeriodicity(periodic, periodic, periodic);
        system.setBoundaries(0.0, sideLength, 0.0, sideLength, 0.0, sideLength);
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

        //        testForce2.addNetwork(hydrogenType, hydrogenType,
        //                              "/home/svenni/Dropbox/projects/programming/fann-md/fann-md/tools/train/tmp/two/fann_network.net",
        //                              "/home/svenni/Dropbox/projects/programming/fann-md/fann-md/tools/train/tmp/two/bounds.fann");


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
        // Minimum for energy conservation: 0.001
        double factor = 1;
        integrator.setTimeStep(factor*1.0);
        system.setIntegrator(&integrator);

        Friction frictionModifier(&system);

        FileManager fileManager(&system);
        fileManager.setOutFileName("/tmp/fannforce/hydrogen-small/atoms*.bin");
        fileManager.setUnitLength(5.2917721092e-11);
        fileManager.setUnitMass(9.10938291e-31);
        system.setFileManager(&fileManager);

        system.setSaveEnabled(true);
        system.setSaveEveryNSteps(1000 / factor);
        system.setOutputEnabled(true);

        system.setupCells(cellCutoff);

        cout << "Cells: " << system.nCells() << endl;


        if(thermo) {
            system.setNSimulationSteps(100000 / factor);
            system.simulate();
//            AndersenThermostat thermostat(&system);
//            thermostat.setCollisionTime(10000);
            BerendsenThermostat thermostat(&system);
            thermostat.setRelaxationTime(100 * factor);
            system.addModifier(&thermostat);
            thermostat.setTargetTemperature(targetTemperature);
            system.setNSimulationSteps(100000 / factor);
            system.simulate();

            system.removeModifier(&thermostat);
            system.setNSimulationSteps(1000000 / factor);
            system.simulate();
        } else if(friction) {
            system.setNSimulationSteps(40000);
            system.simulate();
            system.addModifier(&frictionModifier);
            system.setNSimulationSteps(20000);
            system.simulate();
            system.removeModifier(&frictionModifier);
            system.setNSimulationSteps(40000);
            system.simulate();
        } else {
            system.setNSimulationSteps(10000);
            system.simulate();

        }

        //        CHECK_EQUAL(3, system.nAtomsTotal());

        LOG(INFO) << "FannForceWater test complete";
    }

    //    TEST(FannForceWater)
    //    {
    //        LOG(INFO) << "FannForceWater test started";

    //        // Particle types
    //        AtomType hydrogenType(0);
    //        hydrogenType.setMass(1.0*1837.0);
    //        hydrogenType.setNumber(1);
    //        hydrogenType.setAbbreviation("H");

    //        AtomType oxygenType(1);
    //        oxygenType.setMass(15.9994*1837.0);
    //        oxygenType.setNumber(8);
    //        oxygenType.setAbbreviation("O");

    //        vector<AtomType> particleTypes;
    //        particleTypes.push_back(hydrogenType);
    //        particleTypes.push_back(oxygenType);

    //        vector<Atom *> atoms;

    //        double sideLength = 9.4;

    //        int type = 2;
    //        if(type == 0)  {
    //            int nx = 3;
    //            int ny = 3;
    //            int nz = 3;
    //            double spacingx = sideLength / nx;
    //            double spacingy = sideLength / ny;
    //            double spacingz = sideLength / nz;
    //            int idCounter = 0;
    //            for(int i = 0; i < nz; i++) {
    //                for (int j = 0; j < ny; j++) {
    //                    for (int k = 0; k < nx; ++k) {
    //                        Atom *hydrogenAtom1 = new Atom(hydrogenType);
    //                        hydrogenAtom1->setPosition(Vector3(-0.7575 + spacingx*k, 0.58707 + spacingy*j, spacingz*i));
    //                        //            hydrogenAtom1->setPosition(Vector3(-0.5, 0.58707, i));
    //                        hydrogenAtom1->setID(1 + idCounter);
    //                        Atom *hydrogenAtom2 = new Atom(hydrogenType);
    //                        hydrogenAtom2->setPosition(Vector3(0.7575 + spacingx*k, 0.58707 + spacingy*j, spacingz*i));
    //                        //            hydrogenAtom2->setPosition(Vector3(0.5, 0.58707, i));
    //                        hydrogenAtom2->setID(2 + idCounter);
    //                        Atom *oxygenAtom = new Atom(oxygenType);
    //                        oxygenAtom->setPosition(Vector3(0.0 + spacingx*k, 0.0 + spacingy*j, spacingz*i));
    //                        oxygenAtom->setID(3 + idCounter);
    //                        atoms.push_back(hydrogenAtom1);
    //                        atoms.push_back(hydrogenAtom2);
    //                        atoms.push_back(oxygenAtom);
    //                        idCounter += 3;
    //                    }
    //                }
    //            }
    //        } else if(type == 1) {
    //            Atom *hydrogenAtom1 = new Atom(hydrogenType);
    //            hydrogenAtom1->setPosition(Vector3(-0.7575, 0.58707, 0));
    //            hydrogenAtom1->setID(1);
    //            Atom *hydrogenAtom2 = new Atom(hydrogenType);
    //            hydrogenAtom2->setPosition(Vector3(0.7575, 0.58707, 0));
    //            //            hydrogenAtom2->setPosition(Vector3(0.5, 0.58707, i));
    //            hydrogenAtom2->setID(2);
    //            Atom *oxygenAtom = new Atom(oxygenType);
    //            oxygenAtom->setPosition(Vector3(0.0, 0.0, 0));
    //            oxygenAtom->setID(3);
    //            atoms.push_back(hydrogenAtom1);
    //            atoms.push_back(hydrogenAtom2);
    //            atoms.push_back(oxygenAtom);

    ////            // Move molecule to center
    ////            for(Atom *atom : atoms) {
    ////                atom->setPosition(atom->position() + Vector3(sideLength / 2.0, sideLength / 2.0, sideLength / 2.0));
    ////            }
    //        } else if(type == 2) {
    //            Atom *hydrogenAtom1 = new Atom(hydrogenType);
    //            hydrogenAtom1->setPosition(Vector3(-0.7575, 1-0.58707, 0));
    //            hydrogenAtom1->setID(1);
    //            Atom *hydrogenAtom2 = new Atom(hydrogenType);
    //            hydrogenAtom2->setPosition(Vector3(0.7575, 1-0.58707, 0));
    //            //            hydrogenAtom2->setPosition(Vector3(0.5, 0.58707, i));
    //            hydrogenAtom2->setID(2);
    //            Atom *oxygenAtom = new Atom(oxygenType);
    //            oxygenAtom->setPosition(Vector3(0.0, 1+0.0, 0));
    //            oxygenAtom->setID(3);
    //            double distance = 3.0;
    //            Atom *hydrogenAtom3 = new Atom(hydrogenType);
    //            hydrogenAtom3->setPosition(Vector3(-0.7575, distance+0.58707, 0));
    //            hydrogenAtom3->setID(4);
    //            Atom *hydrogenAtom4 = new Atom(hydrogenType);
    //            hydrogenAtom4->setPosition(Vector3(0.7575, distance+0.58707, 0));
    //            hydrogenAtom4->setID(5);
    //            Atom *oxygenAtom2 = new Atom(oxygenType);
    //            oxygenAtom2->setPosition(Vector3(0.0, distance+0.0, 0));
    //            oxygenAtom2->setID(6);
    //            atoms.push_back(hydrogenAtom1);
    //            atoms.push_back(hydrogenAtom2);
    //            atoms.push_back(oxygenAtom);
    //            atoms.push_back(hydrogenAtom3);
    //            atoms.push_back(hydrogenAtom4);
    //            atoms.push_back(oxygenAtom2);
    //        }

    //        bool thermo = false;
    //        if(thermo) {
    //            Generator gen;
    //            gen.boltzmannDistributeVelocities(0.001, atoms);
    //        }

    //        MoleculeSystem system;
    //        system.setPeriodicity(true, true, true);
    //        system.setBoundaries(0.0, sideLength, 0.0, sideLength, 0.0, sideLength);
    //        system.setParticleTypes(particleTypes);
    //        system.addAtoms(atoms);

    //        FannTwoParticleForce testForce2;
    //        testForce2.addNetwork(hydrogenType, hydrogenType,
    //                              "/home/svenni/Dropbox/studies/master/results/fann_train/20140419-164116/fann_network.net",
    //                              "/home/svenni/Dropbox/studies/master/results/fann_train/20140419-164116/bounds.fann");
    //        testForce2.addNetwork(hydrogenType, oxygenType,
    //                              "/home/svenni/Dropbox/studies/master/results/fann_train/20140418-165017/fann_network.net",
    //                              "/home/svenni/Dropbox/studies/master/results/fann_train/20140418-165017/bounds.fann");
    //        testForce2.addNetwork(oxygenType, oxygenType,
    //                              "/home/svenni/Dropbox/studies/master/results/fann_train/20140418-162313/fann_network.net",
    //                              "/home/svenni/Dropbox/studies/master/results/fann_train/20140418-162313/bounds.fann");

    //        //        LennardJonesForce testForce2;
    //        //        testForce2.setPotentialConstant(pow(2, -1.0/6.0) * 5 / 1.89);
    //        //        testForce2.setShift(3.6 / 1.89);
    //        //        testForce2.setEnergyConstant(1);

    //        testForce2.setCutoffRadius(sideLength / 3);

    //        cout << "Cutoff: " << testForce2.cutoffRadius() << endl;

    //        system.setTwoParticleForce(&testForce2);

    //        FannThreeParticleForce testForce3;
    //        bool flat = false;
    //        if(flat) {
    //            testForce3.loadNetwork("/home/svenni/Dropbox/studies/master/results/fann_train/20140420-220410/fann_network.net",
    //                                   "/home/svenni/Dropbox/studies/master/results/fann_train/20140420-220410/bounds.fann");
    //        } else {
    //            testForce3.loadNetwork("/home/svenni/Dropbox/studies/master/results/fann_train/20140420-001516/fann_network.net",
    //                                   "/home/svenni/Dropbox/studies/master/results/fann_train/20140420-001516/bounds.fann");
    //        }
    //        system.setThreeParticleForce(&testForce3);

    //        VelocityVerletIntegrator integrator(&system);
    //        integrator.setTimeStep(0.01);
    //        system.setIntegrator(&integrator);

    //        BerendsenThermostat thermostat(&system);
    //        if(thermo) {
    //            system.addModifier(&thermostat);
    //        }

    //        FileManager fileManager(&system);
    //        fileManager.setOutFileName("/tmp/fannforce/water/atoms*.xyz");
    //        fileManager.setUnitLength(1.0e-10);
    //        system.setFileManager(&fileManager);

    //        system.setSaveEnabled(true);
    //        system.setSaveEveryNSteps(100);
    //        system.setOutputEnabled(true);

    //        system.setupCells();


    //        if(thermo) {
    //            system.setNSimulationSteps(10000);

    //            thermostat.setTargetTemperature(0.01);
    //            system.simulate();
    //            thermostat.setTargetTemperature(0.005);
    //            system.simulate();
    //            thermostat.setTargetTemperature(0.001);
    //            system.simulate();
    //            thermostat.setTargetTemperature(0.0001);
    //            system.simulate();
    //            system.removeModifier(&thermostat);
    //            system.simulate();
    //        } else {
    //            system.setNSimulationSteps(100000);
    //            system.simulate();
    //        }

    //        //        CHECK_EQUAL(3, system.nAtomsTotal());

    //        LOG(INFO) << "FannForceWater test complete";
    //    }
}
