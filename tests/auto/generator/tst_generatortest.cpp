#include <src/generator.h>
#include <src/atom.h>

#include <QString>
#include <QtTest>

class GeneratorTest : public QObject
{
    Q_OBJECT
    
public:
    GeneratorTest();
    
private Q_SLOTS:
    void initTestCase();
    void cleanupTestCase();
    void boltzmannDistributeVelocities();
};

GeneratorTest::GeneratorTest()
{
}

void GeneratorTest::initTestCase()
{
}

void GeneratorTest::cleanupTestCase()
{
}

/*!
 * \brief GeneratorTest::boltzmannDistributeVelocities ensures that the velocity distribution has a total velocity drift of zero.
 */
void GeneratorTest::boltzmannDistributeVelocities()
{
    Generator generator;
    vector<Atom*> atoms = generator.generateFcc(2,8,AtomType::argon());
    generator.boltzmannDistributeVelocities(100, atoms);

    Vector3 totalVelocity; // = zeros<Vector3>(3);
    for(Atom* atom : atoms) {
        totalVelocity += atom->velocity();
    }
    double maxVelocity = 1;
    QVERIFY(fabs(maxVelocity) < 1e-5);
}

QTEST_APPLESS_MAIN(GeneratorTest)

#include "tst_generatortest.moc"
