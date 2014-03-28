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
    vector<Atom*> atoms = generator.generateFcc(bUnit, 6, AtomType::argon());
//    generator.boltzmannDistributeVelocities(20.0, atoms);

    VelocityVerletIntegrator integrator(&system);
    integrator.setTimeStep(0.005);
    system.setIntegrator(&integrator);
//    system.setInteratomicForce(&force);
    system.setTwoParticleForce(&force);
//    system.setPotentialConstant(potentialConstant);
    system.setBoundaries(generator.lastBoundaries());
    system.addAtoms(atoms);
    system.setupCells();
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
        atom1->addPotential(1.0 / 3.0);
        atom2->addPotential(1.0 / 3.0);
        atom3->addPotential(1.0 / 3.0);
    }
    void calculateAndApplyForce(Atom *atom1, Atom *atom2, Atom *atom3, const Vector3 &atom2Offset, const Vector3 &atom3Offset)
    {
        (void)atom1;
        (void)atom2;
        (void)atom3;
        (void)atom2Offset;
        (void)atom3Offset;
        atom1->addPotential(1.0 / 3.0);
        atom2->addPotential(1.0 / 3.0);
        atom3->addPotential(1.0 / 3.0);
    }
};

TEST(ThreeParticleForceSystem)
{
    cout << "Testing three-particle forces" << endl;
    double potentialConstant = 1;
    MoleculeSystem system;
    system.setOutputEnabled(true);
    system.setSaveEnabled(true);
    Generator generator;
    TwoParticleTestForce testForce2;
    testForce2.setCutoffRadius(1.0);
    ThreeParticleTestForce testForce3;
    vector<Atom*> atoms = generator.generateFcc(1.0, 5, AtomType::argon());
    cout << "N atoms: " << atoms.size() << endl;

    VelocityVerletIntegrator integrator(&system);
    integrator.setTimeStep(0.005);
    system.setIntegrator(&integrator);
    system.setTwoParticleForce(&testForce2);
    system.setThreeParticleForce(&testForce3);
    system.setBoundaries(generator.lastBoundaries());
    system.addAtoms(atoms);
    system.setSaveEnabled(false);
    system.setSaveEveryNSteps(10);
    system.setOutputEnabled(true);
    system.setNSimulationSteps(100);
    system.setupCells();

    testForce2.setCutoffRadius(1000.0);

    system.simulate();
}
