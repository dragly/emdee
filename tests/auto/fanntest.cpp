#include <unittest++/UnitTest++.h>

#include <moleculesystem.h>
#include <moleculesystemcell.h>
#include <generator.h>
#include <atomtype.h>
#include <force/lennardjonesforce.h>
#include <force/fannthreeparticleforce.h>
#include <force/fanntwoparticleforce.h>
#include <modifier/andersenthermostat.h>
#include <atom.h>
#include <processor.h>
#include <utils/logging.h>
#include <force/threeparticleforce.h>
#include <integrator/velocityverletintegrator.h>
#include <filemanager.h>

SUITE(FannForceSystem) {
    TEST(Dummy) {

    }

//    TEST(FannFunTest) {
//        // Particle types
//        AtomType hydrogenType(0);
//        hydrogenType.setMass(1.0);
//        hydrogenType.setNumber(1);
//        hydrogenType.setAbbreviation("H");
//        ofstream outfile("data.out");
//        for(int i = 0; i < 100; i++) {
//            Atom *hydrogenAtom1 = new Atom(hydrogenType);
//            hydrogenAtom1->setPosition(Vector3(0,0,0));
//            Atom *hydrogenAtom2 = new Atom(hydrogenType);
//            hydrogenAtom2->setPosition(Vector3(0,1.0/1.89+ double(i)/100.0 * 6/1.89,0));
//            FannTwoParticleForce testForce2;
//            testForce2.setNewtonsThirdLawEnabled(true);
//            testForce2.loadNetwork("/home/svenni/Dropbox/studies/master/results/fann_train/20140418-165800/fann_network.net",
//                                   "/home/svenni/Dropbox/studies/master/results/fann_train/20140418-165800/bounds.fann");
//            testForce2.setCutoffRadius(6.0 / 1.89);
//            testForce2.calculateAndApplyForce(hydrogenAtom1, hydrogenAtom2);
//            outfile << hydrogenAtom2->position().y()*1.8897 << " " << hydrogenAtom1->force().y() << endl;
//        }
//        outfile.close();
//    }

    TEST(FannForceWater)
    {
        LOG(INFO) << "FannForceWater test started";

        // Particle types
        AtomType hydrogenType(0);
        hydrogenType.setMass(1.0);
        hydrogenType.setNumber(1);
        hydrogenType.setAbbreviation("H");

        AtomType oxygenType(1);
        oxygenType.setMass(15.9994);
        oxygenType.setNumber(8);
        oxygenType.setAbbreviation("O");

        vector<AtomType> particleTypes;
        particleTypes.push_back(hydrogenType);
        particleTypes.push_back(oxygenType);

        vector<Atom *> atoms;
        for(int i = 0; i < 1; i++) {
            Atom *hydrogenAtom1 = new Atom(hydrogenType);
            hydrogenAtom1->setPosition(Vector3(-0.7575, 0.58707, 2*i));
//            hydrogenAtom1->setPosition(Vector3(-0.5, 0.58707, i));
            hydrogenAtom1->setID(1 + i*3);
            Atom *hydrogenAtom2 = new Atom(hydrogenType);
            hydrogenAtom2->setPosition(Vector3(0.7575, 0.58707, 2*i));
//            hydrogenAtom2->setPosition(Vector3(0.5, 0.58707, i));
            hydrogenAtom2->setID(2 + i*3);
            Atom *oxygenAtom = new Atom(oxygenType);
            oxygenAtom->setPosition(Vector3(0.0, 0.0, 2*i));
            oxygenAtom->setID(3 + i*3);
            atoms.push_back(hydrogenAtom1);
            atoms.push_back(hydrogenAtom2);
            atoms.push_back(oxygenAtom);
        }

//        for(int i = 0; i < 5; i++) {
//            for(int j = 0; j < 5; j++) {
//                for(int k = 0; k < 5; k++) {
//                    Atom *hydrogenAtom1 = new Atom(hydrogenType);
//                    hydrogenAtom1->setPosition(Vector3(0,1.0/1.89*j,1.0/1.89*k));
//                    //            hydrogenAtom1->setPosition(Vector3(i*1.0, 0.0, 0.0));
//                    hydrogenAtom1->setID(/*i*5*5 + */j*5 + k + 1);
//                    //            hydrogenAtom1->setPosition(Vector3(0.0, 0.0, 0.0));
//                    //            hydrogenAtom1->setID(1);
//                    //            Atom *hydrogenAtom2 = new Atom(hydrogenType);
//                    //            hydrogenAtom2->setPosition(Vector3(0.74 + 1, i*0.8 + 1, 0.0));
//                    //            hydrogenAtom2->setPosition(Vector3(i*1.0, 0.74, 0.0));
//                    //            hydrogenAtom2->setID(i* 2 + 2);
//                    //            hydrogenAtom2->setPosition(Vector3(0.0, 0.0, 3.0));
//                    //            hydrogenAtom2->setID(2);
//                    //        Atom oxygenAtom(oxygenType);
//                    //        oxygenAtom.setPosition(Vector3(0.0, 0.0, 0.0));
//                    //        oxygenAtom.setID(3);

//                    atoms.push_back(hydrogenAtom1);
//                }
//            }
//        }
        //        atoms.push_back(&oxygenAtom);

        // Move molecule to center
        for(Atom *atom : atoms) {
            atom->setPosition(atom->position() + Vector3(10.0, 10.0, 10.0));
        }

        MoleculeSystem system;
        system.setPeriodicity(true, true, true);
        system.setBoundaries(0.0, 50.0, 0.0, 50.0, 0.0, 50.0);
        system.setParticleTypes(particleTypes);
        system.addAtoms(atoms);

        FannTwoParticleForce testForce2;
        testForce2.addNetwork(hydrogenType, hydrogenType,
                               "/home/svenni/Dropbox/studies/master/results/fann_train/20140419-164116/fann_network.net",
                               "/home/svenni/Dropbox/studies/master/results/fann_train/20140419-164116/bounds.fann");
        testForce2.addNetwork(hydrogenType, oxygenType,
                               "/home/svenni/Dropbox/studies/master/results/fann_train/20140418-165017/fann_network.net",
                               "/home/svenni/Dropbox/studies/master/results/fann_train/20140418-165017/bounds.fann");
        testForce2.addNetwork(oxygenType, oxygenType,
                               "/home/svenni/Dropbox/studies/master/results/fann_train/20140418-162313/fann_network.net",
                               "/home/svenni/Dropbox/studies/master/results/fann_train/20140418-162313/bounds.fann");

//        LennardJonesForce testForce2;
//        testForce2.setPotentialConstant(pow(2, -1.0/6.0) * 5 / 1.89);
//        testForce2.setShift(3.6 / 1.89);
//        testForce2.setEnergyConstant(1);

        testForce2.setCutoffRadius(6.0 / 1.89);
        system.setTwoParticleForce(&testForce2);

        FannThreeParticleForce testForce3;
        testForce3.loadNetwork("/home/svenni/Dropbox/studies/master/results/fann_train/20140420-001516/fann_network.net",
                               "/home/svenni/Dropbox/studies/master/results/fann_train/20140420-001516/bounds.fann");
//        testForce3.loadNetwork("/home/svenni/Dropbox/projects/programming/fann-md/fann-md/tools/train/tmp/fann_network.net",
//                               "/home/svenni/Dropbox/projects/programming/fann-md/fann-md/tools/train/tmp/bounds.fann");
        system.setThreeParticleForce(&testForce3);

        VelocityVerletIntegrator integrator(&system);
        integrator.setTimeStep(0.001);
        system.setIntegrator(&integrator);

        //        AndersenThermostat thermostat(&system);
        //        thermostat.setTargetTemperature(0.01);
        //        system.addModifier(&thermostat);

        FileManager fileManager(&system);
        fileManager.setOutFileName("/tmp/fannforce/atoms*.xyz");
        fileManager.setUnitLength(1.0e-10);
        system.setFileManager(&fileManager);

        system.setSaveEnabled(true);
        system.setSaveEveryNSteps(100);
        system.setOutputEnabled(true);
        system.setNSimulationSteps(100000);
        system.setupCells();
        system.simulate();

        //        CHECK_EQUAL(3, system.nAtomsTotal());

        LOG(INFO) << "FannForceWater test complete";
    }
}
