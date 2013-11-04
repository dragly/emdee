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

# The .cpp file which was generated for your project. Feel free to hack it.
SOURCES += main.cpp \
    moleculardynamics.cpp \
    supercamera.cpp \
    stereoviewport.cpp \
    qglmaskedsurface.cpp \
    qgldrawbuffersurface.cpp \
    fileio.cpp \
    oculusreader.cpp

# Oculus libraries
linux: INCLUDEPATH += /home/svenni/apps/oculus/Oculus/OculusSDK/LibOVR/Include
linux: LIBS += -L/home/svenni/apps/oculus/Oculus/OculusSDK/LibOVR/Lib/Linux/Release/x86_64/ -lovr
linux: LIBS += -lX11
linux: LIBS += -lXinerama
linux: LIBS += -ludev

# Installation path
# target.path =

# Please do not modify the following two lines. Required for deployment.
include(qtquick2applicationviewer/qtquick2applicationviewer.pri)
qtcAddDeployment()

HEADERS += \
    moleculardynamics.h \
    supercamera.h \
    stereoviewport.h \
    qglmaskedsurface_p.h \
    qgldrawbuffersurface_p.h \
    fileio.h \
    oculusreader.h

include(../molecular-dynamics.pri)
message(find $$SRC_DIR -name '*.cpp')
SOURCES += $$system(find $$SRC_DIR -name \'*.cpp\')
SOURCES = $$replace(SOURCES, $$SRC_DIR/main.cpp, )
SOURCES -= $$SRC_DIR/processor.cpp
SOURCES -= $$SRC_DIR/filemanager.cpp
SOURCES -= $$SRC_DIR/configurationparser.cpp

#message($$SOURCES)
INCLUDEPATH += /home/svenni/apps/armadillo

OTHER_FILES += \
    android/src/org/kde/necessitas/ministro/IMinistro.aidl \
    android/src/org/kde/necessitas/ministro/IMinistroCallback.aidl \
    android/src/org/qtproject/qt5/android/bindings/QtActivity.java \
    android/src/org/qtproject/qt5/android/bindings/QtApplication.java \
    android/AndroidManifest.xml \
    android/res/values-fa/strings.xml \
    android/res/values-ru/strings.xml \
    android/res/values-et/strings.xml \
    android/res/values/strings.xml \
    android/res/values/libs.xml \
    android/res/values-es/strings.xml \
    android/res/values-id/strings.xml \
    android/res/values-pl/strings.xml \
    android/res/values-pt-rBR/strings.xml \
    android/res/values-fr/strings.xml \
    android/res/values-nb/strings.xml \
    android/res/layout/splash.xml \
    android/res/values-rs/strings.xml \
    android/res/values-ms/strings.xml \
    android/res/values-ja/strings.xml \
    android/res/values-de/strings.xml \
    android/res/values-zh-rCN/strings.xml \
    android/res/values-nl/strings.xml \
    android/res/values-ro/strings.xml \
    android/res/values-it/strings.xml \
    android/res/values-el/strings.xml \
    android/res/values-zh-rTW/strings.xml \
    android/version.xml \
    android/src/org/kde/necessitas/ministro/IMinistro.aidl \
    android/src/org/kde/necessitas/ministro/IMinistroCallback.aidl \
    android/src/org/qtproject/qt5/android/bindings/QtActivity.java \
    android/src/org/qtproject/qt5/android/bindings/QtApplication.java \
    android/AndroidManifest.xml \
    android/res/values-fa/strings.xml \
    android/res/values-ru/strings.xml \
    android/res/values-et/strings.xml \
    android/res/values/strings.xml \
    android/res/values/libs.xml \
    android/res/values-es/strings.xml \
    android/res/values-id/strings.xml \
    android/res/values-pl/strings.xml \
    android/res/values-pt-rBR/strings.xml \
    android/res/values-fr/strings.xml \
    android/res/values-nb/strings.xml \
    android/res/layout/splash.xml \
    android/res/values-rs/strings.xml \
    android/res/values-ms/strings.xml \
    android/res/values-ja/strings.xml \
    android/res/values-de/strings.xml \
    android/res/values-zh-rCN/strings.xml \
    android/res/values-nl/strings.xml \
    android/res/values-ro/strings.xml \
    android/res/values-it/strings.xml \
    android/res/values-el/strings.xml \
    android/res/values-zh-rTW/strings.xml \
    android/version.xml
