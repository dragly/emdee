#include <unittest++/UnitTest++.h>

#include <moleculesystem.h>
#include <moleculesystemcell.h>
#include <generator.h>
#include <atomtype.h>
#include <force/lennardjonesforce.h>
#include <force/fannforce.h>
#include <force/fanntwoparticleforce.h>
#include <atom.h>
#include <processor.h>
#include <utils/logging.h>
#include <force/threeparticleforce.h>
#include <integrator/velocityverletintegrator.h>
#include <filemanager.h>
#include <force/fannforce.h>
#include <force/fanntwoparticleforce.h>

SUITE(FannForceSystem) {
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

        Atom hydrogenAtom1(hydrogenType);
        hydrogenAtom1.setPosition(Vector3(-0.81411, 0.5807, 0.0));
//        hydrogenAtom1.setPosition(Vector3(-0.4, 0.2, 0.0));
        hydrogenAtom1.setID(1);
        Atom hydrogenAtom2(hydrogenType);
                hydrogenAtom2.setPosition(Vector3(0.81411, 0.5807, 0.0));
//        hydrogenAtom2.setPosition(Vector3(0.5, 0.5, 0.0));
        hydrogenAtom2.setID(2);
        Atom oxygenAtom(oxygenType);
        oxygenAtom.setPosition(Vector3(0.0, 0.0, 0.0));
        oxygenAtom.setID(3);

        vector<Atom *> atoms;
        atoms.push_back(&hydrogenAtom1);
        atoms.push_back(&hydrogenAtom2);
        atoms.push_back(&oxygenAtom);

        // Move molecule to center
        for(Atom *atom : atoms) {
            atom->setPosition(atom->position() + Vector3(5.0, 5.0, 5.0));
        }

        MoleculeSystem system;
        system.setPeriodicity(false, false, false);
        system.setBoundaries(0.0, 10.0, 0.0, 10.0, 0.0, 10.0);
        system.setParticleTypes(particleTypes);
        system.addAtoms(atoms);

        FannTwoParticleForce testForce2;
        testForce2.setCutoffRadius(3.0);

        FannForce testForce3;
        testForce3.loadNetwork("/home/svenni/Dropbox/projects/programming/fann-md/build-fann-md-hyperon-Desktop_Qt_5_2_1_GCC_64bit-Release/fann-trainer/fann_network0.net");
        system.setTwoParticleForce(&testForce2);
        system.setThreeParticleForce(&testForce3);

        VelocityVerletIntegrator integrator(&system);
        integrator.setTimeStep(0.01);
        system.setIntegrator(&integrator);

        FileManager fileManager(&system);
        fileManager.setOutFileName("/tmp/fannforce/atoms*.xyz");
        fileManager.setUnitLength(1.0e-10);
        system.setFileManager(&fileManager);

        system.setSaveEnabled(true);
        system.setSaveEveryNSteps(50);
        system.setOutputEnabled(true);
        system.setNSimulationSteps(10000);
        system.setupCells();
        system.simulate();

        CHECK_EQUAL(3, system.nAtomsTotal());

        LOG(INFO) << "FannForceWater test complete";
    }
}
