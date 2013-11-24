QT += 3dquick
# Add more folders to ship with the application, here
folder_01.source = qml/gui
folder_01.target = qml
DEPLOYMENTFOLDERS = folder_01

# Additional import path used to resolve QML modules in Creator's code model
QML_IMPORT_PATH =

# If your application uses the Qt Mobility libraries, uncomment the following
# lines and add the respective components to the MOBILITY variable.
# CONFIG += mobility
# MOBILITY +=

# Installation path
# target.path =

# Please do not modify the following two lines. Required for deployment.
include(qtquick2applicationviewer/qtquick2applicationviewer.pri)
qtcAddDeployment()

include(../../defaults.pri)

TARGET = emdee

# Local libraries
LIBS += -L../libs -lemdee

# Oculus libraries
oculus {
    linux: INCLUDEPATH += /home/svenni/apps/oculus/Oculus/OculusSDK/LibOVR/Include
    linux: LIBS += -L/home/svenni/apps/oculus/Oculus/OculusSDK/LibOVR/Lib/Linux/Release/x86_64/ -lovr
    linux: LIBS += -lX11
    linux: LIBS += -lXinerama
    linux: LIBS += -ludev
}

INCLUDEPATH += /home/svenni/apps/armadillo

# The .cpp file which was generated for your project. Feel free to hack it.
HEADERS += \
    moleculardynamics.h \
    supercamera.h \
    stereoviewport.h \
    qglmaskedsurface_p.h \
    qgldrawbuffersurface_p.h \
    fileio.h

SOURCES += main.cpp \
    moleculardynamics.cpp \
    supercamera.cpp \
    stereoviewport.cpp \
    qglmaskedsurface.cpp \
    qgldrawbuffersurface.cpp \
    fileio.cpp

oculus {
    DEFINES += MD_USE_OCULUS

    HEADERS +=  \
        oculusreader.h

    SOURCES +=  \
        oculusreader.cpp
}
