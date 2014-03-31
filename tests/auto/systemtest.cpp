#include <moleculesystem.h>
#include <moleculesystemcell.h>
#include <generator.h>
#include <atomtype.h>
#include <force/lennardjonesforce.h>
#include <force/fannforce.h>
#include <force/fanntwoparticleforce.h>
#include <atom.h>

#include <force/threeparticleforce.h>
#include <integrator/velocityverletintegrator.h>
#include <filemanager.h>

#include <unittest++/UnitTest++.h>

TEST(CellSetup)
{
    double potentialConstant = 1.0;
    double bUnit = 5.620 / 3.405;
    MoleculeSystem system;
    system.setOutputEnabled(true);
    system.setSaveEnabled(true);
    system.setPeriodicity(true, true, true);
    Generator generator;
//    generator.setUnitLength(unitLength);
    LennardJonesForce force;
    force.setPotentialConstant(potentialConstant);
    vector<AtomType> particleTypes;
    AtomType dummyType(0);
    dummyType.setNumber(1);
    particleTypes.push_back(dummyType);
    system.setParticleTypes(particleTypes);
    vector<Atom*> atoms = generator.generateFcc(bUnit, 6, dummyType);
    int nAtomsToBeginWith = atoms.size();
//    generator.boltzmannDistributeVelocities(20.0, atoms);

    VelocityVerletIntegrator integrator(&system);
    integrator.setTimeStep(0.001);
    system.setIntegrator(&integrator);
//    system.setInteratomicForce(&force);
    system.setTwoParticleForce(&force);
//    system.setPotentialConstant(potentialConstant);
    system.setBoundaries(generator.lastBoundaries());
    system.addAtoms(atoms);
    system.setupCells();
//    system.setUnitLength(unitLength);

    system.setSaveEnabled(false);
    system.setOutputEnabled(false);
    system.setNSimulationSteps(100);
    system.simulate();
    system.refreshCellContents();
    CHECK_EQUAL(nAtomsToBeginWith, system.nAtomsTotal());
}
