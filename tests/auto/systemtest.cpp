#include <moleculesystem.h>
#include <moleculesystemcell.h>
#include <generator.h>
#include <atomtype.h>
#include <force/lennardjonesforce.h>
#include <integrator/velocityverletintegrator.h>
#include <filemanager.h>

#include <unittest++/UnitTest++.h>

TEST(CellSetup)
{
//    double unitLength = 3.405;
    double potentialConstant = 1;
    double bUnit = 5.620 / 3.405;
    MoleculeSystem system;
    system.setOutputEnabled(true);
    system.setSaveEnabled(true);
    FileManager fileManager(&system);
    fileManager.setOutFileName("/tmp/out*.bin");
    system.setFileManager(&fileManager);
    Generator generator;
//    generator.setUnitLength(unitLength);
    LennardJonesForce force;
    force.setPotentialConstant(potentialConstant);
    vector<Atom*> atoms = generator.generateFcc(bUnit, 12, AtomType::argon());
//    generator.boltzmannDistributeVelocities(20.0, atoms);

    VelocityVerletIntegrator integrator(&system);
    integrator.setTimeStep(0.005);
    system.setIntegrator(&integrator);
//    system.setInteratomicForce(&force);
//    system.setPotentialConstant(potentialConstant);
    system.setBoundaries(generator.lastBoundaries());
    system.addAtoms(atoms);
    system.setupCells(potentialConstant * 3);
//    system.setUnitLength(unitLength);

    ulong nMoleculesInCells = 0;
    for(MoleculeSystemCell* cell : system.cells()) {
        nMoleculesInCells += cell->atoms().size();
    }
    CHECK_EQUAL(system.atoms().size(), nMoleculesInCells);
    system.setSaveEnabled(false);
    system.setOutputEnabled(false);
    system.setNSimulationSteps(100);
    system.simulate();
    system.refreshCellContents();
    nMoleculesInCells = 0;
    for(MoleculeSystemCell* cell : system.cells()) {
        nMoleculesInCells += cell->atoms().size();
    }
    CHECK_EQUAL(system.atoms().size(), nMoleculesInCells);
}
