#include <QString>
#include <QtTest>
#include "moleculesystem.h"

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
    string test = "blah test banana.xkd";

    cout << test.substr(test.length() - 4, 4) << endl;

    QVERIFY2(true, "Failure");
}

QTEST_APPLESS_MAIN(PlaygroundTest)

#include "tst_playgroundtest.moc"
