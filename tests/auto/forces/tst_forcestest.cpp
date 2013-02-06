#include <QString>
#include <QtTest>

class ForcesTest : public QObject
{
    Q_OBJECT
    
public:
    ForcesTest();
    
private Q_SLOTS:
    void testCase1();
};

ForcesTest::ForcesTest()
{
}

void ForcesTest::testCase1()
{
    QVERIFY2(true, "Failure");
}

QTEST_APPLESS_MAIN(ForcesTest)

#include "tst_forcestest.moc"
