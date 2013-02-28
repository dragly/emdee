#include <src/moleculesystem.h>
#include <src/moleculesystemcell.h>
#include <src/generator.h>
#include <src/atomtype.h>
#include <src/interatomicforce.h>
#include <src/integrator/velocityverletintegrator.h>

#include <QString>
#include <QtTest>

class SystemTest : public QObject
{
    Q_OBJECT
    
public:
    SystemTest();
    
private Q_SLOTS:
    void initTestCase();
    void cleanupTestCase();
    void cellSetup();
};

SystemTest::SystemTest()
{
}

void SystemTest::initTestCase()
{
}

void SystemTest::cleanupTestCase()
{
}

void SystemTest::cellSetup()
{
    double unitLength = 3.405;
    double potentialConstant = 1;
    double bUnit = 5.620 / 3.405;
    MoleculeSystem system;
    system.setSaveEnabled(false);
    Generator generator;
//    generator.setUnitLength(unitLength);
    InteratomicForce force;
    force.setPotentialConstant(potentialConstant);
    vector<Atom*> atoms = generator.generateFcc(bUnit, 7, AtomType::argon());
    generator.boltzmannDistributeVelocities(20.0, atoms);

    VelocityVerletIntegrator integrator(&system);
    integrator.setTimeStep(0.005);
    system.setIntegrator(&integrator);
    system.setInteratomicForce(&force);
//    system.setPotentialConstant(potentialConstant);
    system.setBoundaries(generator.lastBoundaries());
    system.addAtoms(atoms);
    system.setupCells(potentialConstant * 3);
//    system.setUnitLength(unitLength);

    ulong nMoleculesInCells = 0;
    for(MoleculeSystemCell* cell : system.cells()) {
        nMoleculesInCells += cell->atoms().size();
    }
    QCOMPARE(system.atoms().size(), nMoleculesInCells);

    system.simulate(50);
    system.refreshCellContents();
    nMoleculesInCells = 0;
    for(MoleculeSystemCell* cell : system.cells()) {
        nMoleculesInCells += cell->atoms().size();
    }
    QCOMPARE(system.atoms().size(), nMoleculesInCells);
}

QTEST_APPLESS_MAIN(SystemTest)

#include "tst_systemtest.moc"
