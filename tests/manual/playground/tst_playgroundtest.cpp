#include <QString>
#include <QtTest>
#include <cmath>
#include <armadillo>
#include "moleculesystem.h"
#include "generator.h"

#include <iostream>

using namespace std;
using namespace arma;

class PlaygroundTest : public QObject
{
    Q_OBJECT
    
public:
    PlaygroundTest();
    int n;

private Q_SLOTS:
    void initTestCase();
    void cleanupTestCase();
    void testCase1();
    void setupCells();
//    void benchmarkMod();
//    void benchmarkIf();
};

PlaygroundTest::PlaygroundTest()
{
    n = 1e7;
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
    Generator generator;
    vector<Molecule*> molecules = generator.generateFcc(nCells, b, AtomType::argon());
    generator.boltzmannDistributeVelocities(molecules);
    MoleculeSystem system;
    system.addMolecules(molecules);
    system.setBoundaries(0, b*nCells);
    system.setupCells(4);
}

//void PlaygroundTest::benchmarkMod() {
//    QBENCHMARK {
//        for(int i = 0; i < n; i++) {
//            vec a = 10 * randn<vec>(3);
//            for(int dim = 0; dim < 3; dim++) {
//                a(dim) = fmod(a(dim) + 4, 4);
//            }
//        }
//    }
//}

//void  PlaygroundTest::benchmarkIf() {


//    QBENCHMARK {
//        for(int i = 0; i < n; i++) {
//            vec a = 10 * randn<vec>(3);
//            for(int dim = 0; dim < 3; dim++) {
//                if(a(dim) > 4) {
//                    a(dim) -= 4;
//                } else if(a(dim) < 0) {
//                    a(dim) += 4;
//                }
//            }
//        }
//    }
//}

QTEST_APPLESS_MAIN(PlaygroundTest)

#include "tst_playgroundtest.moc"
