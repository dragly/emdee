#include <src/force/lennardjonesforce.h>
#include <src/atom.h>

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
    Atom* atom1 = new Atom(AtomType::argon());
    Atom* atom2 = new Atom(AtomType::argon());
    rowvec position1;
    position1 << -1 << 0 << 0;
    rowvec position2;
    position2 << 1 << 0 << 0;
    atom1->setPosition(position1);
    atom2->setPosition(position2);
    LennardJonesForce force;
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
