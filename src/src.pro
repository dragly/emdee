include(../molecular-dynamics.pri)

TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

INCLUDEPATH += ../

TARGET = molecular-dynamics

SOURCES += main.cpp \
    moleculesystem.cpp \
    atom.cpp \
    generator.cpp \
    atomtype.cpp \
    molecule.cpp \
    moleculesystemcell.cpp \
    interatomicforce.cpp \
    integrator/velocityverletintegrator.cpp \
    integrator/integrator.cpp \
    integrator/eulercromerintegrator.cpp \
    filemanager.cpp

HEADERS += \
    moleculesystem.h \
    atom.h \
    generator.h \
    atomtype.h \
    molecule.h \
    moleculesystemcell.h \
    interatomicforce.h \
    integrator/velocityverletintegrator.h \
    integrator/integrator.h \
    integrator/eulercromerintegrator.h \
    filemanager.h

# Building
myscript.target = myscript
myscript.commands = python $$PWD/../myscript.py $$PWD/../
QMAKE_EXTRA_TARGETS += myscript
PRE_TARGETDEPS += myscript

OTHER_FILES += ../testconfig.cfg \
    ../molecular-dynamics.pri
