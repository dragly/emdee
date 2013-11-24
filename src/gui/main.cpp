#include <QtGui/QGuiApplication>
#include <QtQuick>
#include "qtquick2applicationviewer.h"
#include <moleculardynamics.h>
#include <stereoviewport.h>
#include <fileio.h>
#include <fenv.h>
#ifdef MD_USE_OCULUS
#include <OVR.h>
#include <oculusreader.h>
#endif

int main(int argc, char *argv[])
{
//    feenableexcept(FE_INVALID | FE_OVERFLOW);
#ifdef MD_USE_OCULUS
    OVR::System::Init();
    qmlRegisterType<OculusReader>("OculusReader", 1, 0, "OculusReader");
#endif
    qmlRegisterType<FileIO>("FileIO", 1, 0, "FileIO");
    qmlRegisterType<StereoViewport>("StereoViewport", 1, 0, "StereoViewport");
    qmlRegisterType<MolecularDynamics>("MolecularDynamics", 1, 0, "MolecularDynamics");

    QGuiApplication app(argc, argv);

    QtQuick2ApplicationViewer viewer;
    viewer.setMainQmlFile(QStringLiteral("qml/gui/main.qml"));
    viewer.showExpanded();

    return app.exec();
}
