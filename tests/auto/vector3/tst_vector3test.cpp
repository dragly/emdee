#include <src/math/vector3.h>

#include <QString>
#include <QtTest>

#include <iostream>

Q_DECLARE_METATYPE(Vector3)

using namespace std;

class Vector3Test : public QObject
{
    Q_OBJECT
    
public:
    Vector3Test();
    
private Q_SLOTS:
    void initTestCase();
    void cleanupTestCase();
    void additionTest();
    void additionTest_data();
    void dotProductTest();
    void dotProductTest_data();
};

Vector3Test::Vector3Test()
{
}

void Vector3Test::initTestCase()
{
}

void Vector3Test::cleanupTestCase()
{
}

void Vector3Test::additionTest()
{
    QFETCH(Vector3, vector1);
    QFETCH(Vector3, vector2);
    QFETCH(Vector3, result);
    Vector3 resulta = vector1 + vector2;

    QCOMPARE(resulta[0], result[0]);
    QCOMPARE(resulta[1], result[1]);
    QCOMPARE(resulta[2], result[2]);
}

void Vector3Test::additionTest_data()
{
    QTest::addColumn<Vector3>("vector1");
    QTest::addColumn<Vector3>("vector2");
    QTest::addColumn<Vector3>("result");
    QTest::newRow("0") << Vector3(1., 2., 3.)
                       << Vector3(4., 5., 6.)
                       << Vector3(5., 7., 9.);
    QTest::newRow("1") << Vector3(1.5, -2.6, 3.1)
                       << Vector3(4., 5.1, -6.)
                       << Vector3(5.5, 2.5, -2.9);
    QTest::newRow("2") << Vector3(1.5, 2.6, 3.1)
                       << Vector3(4., 5.1, 6.)
                       << Vector3(5.5, 7.7, 9.1);
    QTest::newRow("2") << Vector3(1., 2., 3.)
                       << Vector3(4., 5., 7.)
                       << Vector3(5., 7., 10.);
}

void Vector3Test::dotProductTest()
{
    QFETCH(Vector3, vector1);
    QFETCH(Vector3, vector2);
    QFETCH(double, result);
    double resulta = vector1 * vector2;

    QCOMPARE(resulta, result);
}

void Vector3Test::dotProductTest_data()
{
    QTest::addColumn<Vector3>("vector1");
    QTest::addColumn<Vector3>("vector2");
    QTest::addColumn<double>("result");
    QTest::newRow("0") << Vector3(1., 2., 3.)
                       << Vector3(4., 5., 6.)
                       << (4. + 10. + 18);
    QTest::newRow("1") << Vector3(1.5, -2.6, 3.1)
                       << Vector3(4., 5.1, -6.)
                       << (1.5 * 4 - 2.6 * 5.1 - 3.1 * 6);
    QTest::newRow("2") << Vector3(1.5, 2.6, 3.1)
                       << Vector3(4., 5.1, 6.)
                       << (1.5 * 4 + 2.6 * 5.1 + 3.1 * 6);
    QTest::newRow("2") << Vector3(1., 2., 3.)
                       << Vector3(4., 5., 7.)
                       << (4. + 10. + 21.);
}

QTEST_APPLESS_MAIN(Vector3Test)

#include "tst_vector3test.moc"
