#include <src/moleculesystem.h>
#include <src/moleculesystemcell.h>
#include <src/generator.h>
#include <src/atomtype.h>
#include <src/force/lennardjonesforce.h>
#include <src/integrator/velocityverletintegrator.h>

#include <unittest++/UnitTest++.h>

TEST(CellSetup)
{
//    double unitLength = 3.405;
    double potentialConstant = 1;
    double bUnit = 5.620 / 3.405;
    MoleculeSystem system;
    system.setSaveEnabled(false);
    Generator generator;
//    generator.setUnitLength(unitLength);
    LennardJonesForce force;
    force.setPotentialConstant(potentialConstant);
    vector<Atom*> atoms = generator.generateFcc(bUnit, 7, AtomType::argon());
    generator.boltzmannDistributeVelocities(20.0, atoms);

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

    system.simulate();
    system.refreshCellContents();
    nMoleculesInCells = 0;
    for(MoleculeSystemCell* cell : system.cells()) {
        nMoleculesInCells += cell->atoms().size();
    }
    CHECK_EQUAL(system.atoms().size(), nMoleculesInCells);
}
