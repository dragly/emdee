#include <src/interatomicforce.h>
#include <src/atom.h>
#include <src/molecule.h>

#include <QString>
#include <QtTest>
#include <armadillo>

using namespace arma;

class ForcesTest : public QObject
{
    Q_OBJECT
    
public:
    ForcesTest();
    
private Q_SLOTS:
    void testForceBetweenArgons();
};

ForcesTest::ForcesTest()
{
}

void ForcesTest::testForceBetweenArgons()
{
    Molecule* molecule1 = new Molecule();
    Molecule* molecule2 = new Molecule();
    Atom* atom1 = new Atom(molecule1, AtomType::argon());
    Atom* atom2 = new Atom(molecule2, AtomType::argon());
    rowvec position1;
    position1 << -1 << 0 << 0;
    rowvec position2;
    position2 << 1 << 0 << 0;
    molecule1->setPosition(position1);
    molecule2->setPosition(position2);
    InteratomicForce force;
    force.setPotentialConstant(3);
    force.calculateAndApplyForce(atom1, atom2);
//    rowvec calculatedForce = force.force();
//    cout << calculatedForce << endl;
//    QVERIFY(true);
//    double val = calculatedForce(0);
//    QCOMPARE(val, 4.0);
}

QTEST_APPLESS_MAIN(ForcesTest)

#include "tst_forcestest.moc"
