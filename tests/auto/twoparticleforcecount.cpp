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

SUITE(TwoParticleForceCount) {

    class TwoParticleTestForce : public TwoParticleForce {
    public:
        virtual void calculateAndApplyForce(Atom *atom1, Atom *atom2)
        {
            calculateAndApplyForce(atom1, atom2, m_zeroVector);
        }
        virtual void calculateAndApplyForce(Atom *atom1, Atom *atom2, const Vector3 &atom2Offset)
        {
            (void)atom2Offset;
            atom1->addPotential(1.0);
            if(isNewtonsThirdLawEnabled()) {
                atom2->addPotential(1.0);
            }
        }
    };

    TEST(SingleAtomSingleCellPeriodic)
    {
        LOG(INFO) << "SingleAtomSingleCellPeriodic test started";
        MoleculeSystem system;
        system.setOutputEnabled(true);
        system.setSaveEnabled(true);
        Generator generator;
        TwoParticleTestForce testForce2;
        testForce2.setCutoffRadius(1000.0); // Let all atoms interact with each other
        vector<AtomType> particleTypes;
        AtomType dummyType(0);
        particleTypes.push_back(dummyType);
        vector<Atom*> atoms;
        Atom* atom = new Atom(dummyType);
        atom->setPosition(Vector3(4.5, 4.5, 4.5));
        atoms.push_back(atom);

        VelocityVerletIntegrator integrator(&system);
        integrator.setTimeStep(0.005);
        system.setPeriodicity(true, true, true);
        system.setParticleTypes(particleTypes);
        system.setIntegrator(&integrator);
        system.setTwoParticleForce(&testForce2);
        system.setBoundaries(0.0, 10.0, 0.0, 10.0, 0.0, 10.0);
        system.addAtoms(atoms);

        CHECK_EQUAL(1, system.atoms().size());

        system.setSaveEnabled(false);
        system.setSaveEveryNSteps(10);
        system.setOutputEnabled(false);
        system.setNSimulationSteps(1);
        system.setupCells(10.0);

        system.simulate();

        // There is 1 atom
        // The same atom is duplicated to 26 neighbor cells, so the total number
        // of potential points gained is
        //   1 * 26 = 26
        CHECK_CLOSE(26, system.potentialEnergyTotal(), 1e-9);
        LOG(INFO) << "SingleAtomSingleCellPeriodic test complete";
    }

    TEST(ManyAtomsSingleCell)
    {
        LOG(INFO) << "ManyAtomsSingleCell test started";
        MoleculeSystem system;
        system.setOutputEnabled(true);
        system.setSaveEnabled(true);
        Generator generator;
        TwoParticleTestForce testForce2;
        testForce2.setCutoffRadius(1000.0); // Let all atoms interact with each other
        vector<AtomType> particleTypes;
        AtomType dummyType(0);
        particleTypes.push_back(dummyType);
        vector<Atom*> atoms = generator.generateFcc(1.0, 2, dummyType);

        // 2*2*2 * 4 = 32 atoms
        CHECK_EQUAL(32, atoms.size());

        VelocityVerletIntegrator integrator(&system);
        integrator.setTimeStep(0.005);
        system.setPeriodicity(false, false, false);
        system.setParticleTypes(particleTypes);
        system.setIntegrator(&integrator);
        system.setTwoParticleForce(&testForce2);
        system.setBoundaries(0.0, 10.0, 0.0, 10.0, 0.0, 10.0);
        system.addAtoms(atoms);
        system.setSaveEnabled(false);
        system.setSaveEveryNSteps(10);
        system.setOutputEnabled(false);
        system.setNSimulationSteps(1);
        system.setupCells(10.0);

        system.simulate();

        // There are 32 atoms
        // Each atom has every other atom as a neighbor, i.e, 31 neighbors, in this cell
        // This results in N(N-1)/2 force calculations, each with 2 potential points
        //   2 * 32 * (32 - 1) / 2 = 992
        CHECK_CLOSE(992, system.potentialEnergyTotal(), 1e-9);
        LOG(INFO) << "ManyAtomsSingleCell test complete";
    }

    TEST(ManyAtomsSingleCellPeriodic)
    {
        LOG(INFO) << "ManyAtomsSingleCellPeriodic test started";
        MoleculeSystem system;
        system.setOutputEnabled(true);
        system.setSaveEnabled(true);
        Generator generator;
        TwoParticleTestForce testForce2;
        vector<AtomType> particleTypes;
        AtomType dummyType(0);
        dummyType.setNumber(1);
        particleTypes.push_back(dummyType);
        vector<Atom*> atoms = generator.generateFcc(1.0, 2, dummyType);

        // 2*2*2 * 4 = 32 atoms
        CHECK_EQUAL(32, atoms.size());

        VelocityVerletIntegrator integrator(&system);
        integrator.setTimeStep(0.0);
        system.setPeriodicity(true, true, true);
        system.setParticleTypes(particleTypes);
        system.setIntegrator(&integrator);
        system.setTwoParticleForce(&testForce2);
        system.setBoundaries(0.0, 2.0, 0.0, 2.0, 0.0, 2.0);
        system.addAtoms(atoms);
        system.setSaveEnabled(false);
        system.setSaveEveryNSteps(10);
        system.setOutputEnabled(false);
        system.setNSimulationSteps(1);
        system.setupCells(10.0);

        // Let all atoms interact with each other
        testForce2.setCutoffRadius(1000.0);
        // There are N=32 atoms
        // Each atom has every other atom as a neighbor, i.e, 31 neighbors, in this cell
        // This results in N(N-1) / 2 force calculations, where both atoms gain 1 potential, i.e. 2 total
        // In addition, there are 26 neighbor cells with 32 atoms
        // Resulting in an additional N*N*26 force calculations where our atom gains 1 potential
        // The total is then
        //   2*32*(32-1)/2 + 1*32*32*26 = 27616
        system.updateForces();
        system.updateStatistics();
        CHECK_CLOSE(27616, system.potentialEnergyTotal(), 1e-9);

        // Test only forces with closest neighbors, should be 12 for FCC
        // Each neighbor gives each atom 1 potential point
        // We have 32 atoms, so we have a total of
        //   1 * 32 * 12 = 384
        testForce2.setCutoffRadius(sqrt(0.5*0.5 + 0.5*0.5) + 0.05);
        system.updateForces();
        for(Atom *atom : system.atoms()) {
            CHECK_CLOSE(12, atom->potential(), 1e-9);
        }
        system.updateStatistics();
        CHECK_CLOSE(384, system.potentialEnergyTotal(), 1e-9);

        LOG(INFO) << "ManyAtomsSingleCellPeriodic test complete";
    }
}

