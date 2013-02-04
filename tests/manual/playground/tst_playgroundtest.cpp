#include <QString>
#include <QtTest>
#include "moleculesystem.h"
#include "generator.h"

#include <iostream>

using namespace std;

class PlaygroundTest : public QObject
{
    Q_OBJECT
    
public:
    PlaygroundTest();

private Q_SLOTS:
    void initTestCase();
    void cleanupTestCase();
    void testCase1();
    void setupCells();
};

PlaygroundTest::PlaygroundTest()
{
}

void PlaygroundTest::initTestCase()
{
}

void PlaygroundTest::cleanupTestCase()
{
}

void PlaygroundTest::testCase1()
{
    MoleculeSystem* system = new MoleculeSystem();
    mat test = zeros(2,3);
    test(1,0) = 10;
    test(1,1) = 10;
    test(1,2) = 10;
    system->setBoundaries(test);
}

void PlaygroundTest::setupCells() {
    int nCells = 5;
    double b = 5.620;
    vector<Molecule*> molecules = Generator::generateFcc(nCells, b, AtomType::argon());
    Generator::boltzmannDistributeVelocities(molecules);
    MoleculeSystem system;
    system.addMolecules(molecules);
    system.setBoundaries(0, b*nCells);
    system.setupCells(4);
}

QTEST_APPLESS_MAIN(PlaygroundTest)

#include "tst_playgroundtest.moc"
