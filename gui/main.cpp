#include <QtGui/QGuiApplication>
#include <QtQuick>
#include "qtquick2applicationviewer.h"
#include <moleculardynamics.h>
#include <stereoviewport.h>
#include <fileio.h>
#include <OVR.h>
#include <oculusreader.h>

int main(int argc, char *argv[])
{
    OVR::System::Init();
    qmlRegisterType<FileIO>("FileIO", 1, 0, "FileIO");
    qmlRegisterType<StereoViewport>("StereoViewport", 1, 0, "StereoViewport");
    qmlRegisterType<OculusReader>("OculusReader", 1, 0, "OculusReader");
    qmlRegisterType<MolecularDynamics>("MolecularDynamics", 1, 0, "MolecularDynamics");

    QGuiApplication app(argc, argv);

    QtQuick2ApplicationViewer viewer;
    viewer.setMainQmlFile(QStringLiteral("qml/gui/main.qml"));
    viewer.showExpanded();

    return app.exec();
}
