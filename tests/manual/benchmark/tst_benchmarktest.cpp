#include "moleculesystem.h"
#include "generator.h"

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
    double temperature = 10.0;
    vector<Molecule*> molecules = generator.generateFcc(bUnit, nCells, AtomType::argon());
    generator.boltzmannDistributeVelocities(temperature, molecules);
    // Set up molecule system
    MoleculeSystem system;
    system.addMolecules(molecules);
    cout << "addded" << endl;
    system.setBoundaries(generator.lastBoundaries());
    cout << "setbounds" << endl;
    system.setupCells();
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
    QTest::newRow("0") << 7 << 20;
    QTest::newRow("1") << 8 << 20;
}

QTEST_APPLESS_MAIN(BenchmarkTest)

#include "tst_benchmarktest.moc"
