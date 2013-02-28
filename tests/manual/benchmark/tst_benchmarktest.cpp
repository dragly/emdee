#include <src/moleculesystem.h>
#include <src/generator.h>
#include <src/force/lennardjonesforce.h>
#include <src/integrator/velocityverletintegrator.h>

#include <QString>
#include <QtTest>
#include <iostream>

using namespace std;

class BenchmarkTest : public QObject
{
    Q_OBJECT
    
public:
    BenchmarkTest();
    
private Q_SLOTS:
    void initTestCase();
    void cleanupTestCase();
    void benchmarkDifferentSizes();
    void benchmarkDifferentSizes_data();
};

BenchmarkTest::BenchmarkTest()
{
}

void BenchmarkTest::initTestCase()
{
}

void BenchmarkTest::cleanupTestCase()
{
}

void BenchmarkTest::benchmarkDifferentSizes()
{
    QFETCH(int, nCells);
    QFETCH(int, nSimulationSteps);
    Generator generator;

    // Generator specific config
    double bUnit = 5.620 / 3.405;
    double temperature = 1.0;
    double potentialConstant = 1.0;
    vector<Atom*> atoms = generator.generateFcc(bUnit, nCells, AtomType::argon());
    generator.boltzmannDistributeVelocities(temperature, atoms);
    // Set up force
    LennardJonesForce interatomicForce;
    interatomicForce.setPotentialConstant(potentialConstant);
    // Set up molecule system
    MoleculeSystem system;
    system.setSaveEnabled(false);
//    system.setOutFileName("/tmp/mdout/test*.bin");
//    system.setOutputEnabled(false);
    // Set up integrator
    VelocityVerletIntegrator integrator(&system);
    integrator.setTimeStep(0.001);
    system.setIntegrator(&integrator);
    // Add molecules
    system.addAtoms(atoms);
    cout << "addded" << endl;
    system.setBoundaries(generator.lastBoundaries());
    cout << "setbounds" << endl;
    system.setupCells(potentialConstant * 3);
    cout << "Setup cells" << endl;
    QBENCHMARK_ONCE {
        system.simulate(nSimulationSteps);
    }
    QVERIFY2(true, "Failure");
}

void BenchmarkTest::benchmarkDifferentSizes_data()
{
    QTest::addColumn<int>("nCells");
    QTest::addColumn<int>("nSimulationSteps");
//    QTest::newRow("0") << 7 << 20;
    QTest::newRow("0") << 8 << 100;
}

QTEST_APPLESS_MAIN(BenchmarkTest)

#include "tst_benchmarktest.moc"
