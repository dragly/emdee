#include <unittest++/UnitTest++.h>

#include <moleculesystem.h>
#include <moleculesystemcell.h>
#include <generator.h>
#include <atomtype.h>
#include <force/lennardjonesforce.h>
#include <force/fannthreeparticleforce.h>
#include <force/fanntwoparticleforce.h>
#include <atom.h>
#include <processor.h>
#include <utils/logging.h>
#include <force/threeparticleforce.h>
#include <integrator/velocityverletintegrator.h>
#include <filemanager.h>

SUITE(ThreeParticleForceSystem) {

    class TwoParticleTestForce : public TwoParticleForce {
    public:
        virtual void calculateAndApplyForce(Atom *atom1, Atom *atom2)
        {
            (void)atom1;
            (void)atom2;
        }
        virtual void calculateAndApplyForce(Atom *atom1, Atom *atom2, const Vector3 &atom2Offset)
        {
            (void)atom1;
            (void)atom2;
            (void)atom2Offset;
        }
    };

    class ThreeParticleTestForce : public ThreeParticleForce {
    public:
        void calculateAndApplyForce(Atom *atom1, Atom *atom2, Atom *atom3)
        {
            (void)atom1;
            (void)atom2;
            (void)atom3;

            double potentialFactor = 1.0 / 3.0;

            atom1->addPotential(potentialFactor);
            atom2->addPotential(potentialFactor);
            atom3->addPotential(potentialFactor);
        }
        void calculateAndApplyForce(Atom *atom1, Atom *atom2, Atom *atom3, const Vector3 &atom2Offset, const Vector3 &atom3Offset)
        {
            (void)atom1;
            (void)atom2;
            (void)atom3;
            (void)atom2Offset;
            (void)atom3Offset;

            double potentialFactor = 1.0 / 3.0;

            atom1->addPotential(potentialFactor);
            atom2->addPotential(potentialFactor);
            atom3->addPotential(potentialFactor);
        }
    };

    TEST(ThreeParticlePotentialSumSingleCell)
    {
        LOG(INFO) << "ThreeParticleForceSum test started";
        MoleculeSystem system;
        system.setOutputEnabled(true);
        system.setSaveEnabled(true);
        Generator generator;
        TwoParticleTestForce testForce2;
        testForce2.setCutoffRadius(1000.0); // Let all atoms interact with each other
        ThreeParticleTestForce testForce3;
        vector<AtomType> particleTypes;
        AtomType dummyType(0);
        particleTypes.push_back(dummyType);
        vector<Atom*> atoms = generator.generateFcc(1.0, 2, dummyType);

        CHECK_EQUAL(32, atoms.size());

        VelocityVerletIntegrator integrator(&system);
        integrator.setTimeStep(0.005);
        system.setPeriodicity(false, false, false);
        system.setParticleTypes(particleTypes);
        system.setIntegrator(&integrator);
        system.setTwoParticleForce(&testForce2);
        system.setThreeParticleForce(&testForce3);
        system.setBoundaries(0.0, 7.5, 0.0, 7.5, 0.0, 7.5);
        system.addAtoms(atoms);
        system.setSaveEnabled(false);
        system.setSaveEveryNSteps(10);
        system.setOutputEnabled(false);
        system.setNSimulationSteps(1);
        system.setupCells(1.0);

        system.simulate();

        // There are 32 atoms
        // Each atom has every other atom as a neighbor, i.e, 31 neighbors
        // There are 31c2 = 465 ways to pick two neighbors for a three-particle force calculation
        // Resulting in 465 * 32 = 14880 three-particle force calculations
        CHECK_CLOSE(14880, system.potentialEnergyTotal(), 1e-9);
        LOG(INFO) << "ThreeParticleForceSum test complete";
    }

    TEST(ThreeParticlePotentialSumPeriodic)
    {
        LOG(INFO) << "ThreeParticleForceSumPeriodic test started";
        MoleculeSystem system;
        Generator generator;
        TwoParticleTestForce testForce2;
        testForce2.setCutoffRadius(1.1);
        ThreeParticleTestForce testForce3;

        vector<AtomType> particleTypes;
        AtomType dummyType(0);
        dummyType.setNumber(1);
        particleTypes.push_back(dummyType);
        vector<Atom*> atoms;
        atoms.reserve(27);
        int idCounter = 1;
        // Set up a system with one atom in each cell
        for(int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; ++j) {
                for (int k = 0; k < 3; ++k) {
                    Atom* atom = new Atom(dummyType);
                    atom->setID(idCounter);
                    atom->setPosition(Vector3(0.5 + k, 0.5 + j, 0.5 + i));
                    atoms.push_back(atom);
                    idCounter += 1;
                }
            }
        }
        system.addAtoms(atoms);

        VelocityVerletIntegrator integrator(&system);
        integrator.setTimeStep(0.005);
        system.setPeriodicity(true, true, true);
        system.setIntegrator(&integrator);
        system.setTwoParticleForce(&testForce2);
        system.setThreeParticleForce(&testForce3);
        system.setBoundaries(0.0, 3.1, 0.0, 3.1, 0.0, 3.1);
        system.setSaveEnabled(true);
        // TODO Add a check to see if the particle types are actually set
        system.setParticleTypes(particleTypes);
        system.setSaveEveryNSteps(1);
        system.setOutputEnabled(false);
        system.setNSimulationSteps(2);
        system.setupCells(1.01);
        system.simulate();

        CHECK_EQUAL(27, system.nAtomsTotal());

        for(MoleculeSystemCell* cell : system.localCells()) {
            for(Atom* atom : cell->atoms()) {
                CHECK_EQUAL(6, atom->neighborAtoms().size());
            }
        }

        for(MoleculeSystemCell* cell : system.localCells()) {
            for(Atom* atom : cell->atoms()) {
                CHECK_CLOSE(15, atom->potential(), 1e-9);
            }
        }
        // There are 27 atoms
        // Each atom has 6 neighbors
        // There are 6c2 = 15 ways to pick two neighbors for a three-particle force calculation
        // Resulting in 15 * 27 = 405 force calculations in total
        CHECK_CLOSE(405, system.potentialEnergyTotal(), 1e-9);
        LOG(INFO) << "ThreeParticleForceSumPeriodic test complete";
    }

    TEST(ThreeParticlePotentialSumPeriodicLarge)
    {
        LOG(INFO) << "ThreeParticleForceSumPeriodic test started";
        MoleculeSystem system;
        Generator generator;
        TwoParticleTestForce testForce2;
        testForce2.setCutoffRadius(1.0);

        ThreeParticleTestForce testForce3;

        vector<AtomType> particleTypes;
        AtomType dummyType(0);
        dummyType.setNumber(1);
        particleTypes.push_back(dummyType);
        vector<Atom*> atoms;
        atoms.reserve(10*10*10);
        int idCounter = 1;
        // Set up a system with one atom in each cell
        for(int i = 0; i < 10; i++) {
            for (int j = 0; j < 10; ++j) {
                for (int k = 0; k < 10; ++k) {
                    Atom* atom = new Atom(dummyType);
                    atom->setID(idCounter);
                    atom->setPosition(Vector3(0.5 + k, 0.5 + j, 0.5 + i));
                    atoms.push_back(atom);
                    idCounter += 1;
                }
            }
        }
        system.addAtoms(atoms);

        VelocityVerletIntegrator integrator(&system);
        integrator.setTimeStep(1.0);
        system.setPeriodicity(true, true, true);
        system.setIntegrator(&integrator);
        system.setTwoParticleForce(&testForce2);
        system.setThreeParticleForce(&testForce3);
        system.setBoundaries(0.0, 10.0, 0.0, 10.0, 0.0, 10.0);
        system.setSaveEnabled(true);
        // TODO Add a check to see if the particle types are actually set
        system.setParticleTypes(particleTypes);
        system.setSaveEveryNSteps(1);
        system.setOutputEnabled(false);
        system.setNSimulationSteps(2);
        system.setupCells(2.0);
        system.simulate();

        CHECK_EQUAL(1000, system.nAtomsTotal());

        for(MoleculeSystemCell* cell : system.localCells()) {
            for(Atom* atom : cell->atoms()) {
                CHECK_EQUAL(6, atom->neighborAtoms().size());
                CHECK_CLOSE(15, atom->potential(), 1e-9);
            }
        }

        // There are 1000 atoms
        // Each atom has 6 neighbors
        // There are 6c2 = 15 ways to pick two neighbors for a three-particle force calculation
        // Resulting in 15 * 1000 = 15000 force calculations in total
        CHECK_CLOSE(15000, system.potentialEnergyTotal(), 1e-9);

        // Enable motion for a couple of time steps
        for(MoleculeSystemCell* cell : system.localCells()) {
            for(Atom* atom : cell->atoms()) {
                atom->setVelocity(Vector3(0.1,0.2,0.3));
            }
        }
        system.setNSimulationSteps(88);
        system.simulate();

        CHECK_EQUAL(1000, system.nAtomsTotal());

        for(MoleculeSystemCell* cell : system.localCells()) {
            for(Atom* atom : cell->atoms()) {
                CHECK_EQUAL(6, atom->neighborAtoms().size());
                CHECK_CLOSE(15, atom->potential(), 1e-9);
            }
        }

        LOG(INFO) << "ThreeParticleForceSumPeriodic test complete";
    }
}
