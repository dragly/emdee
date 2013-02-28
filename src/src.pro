include(../molecular-dynamics.pri)

TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

INCLUDEPATH += ../

TARGET = molecular-dynamics

SOURCES += main.cpp \
    moleculesystem.cpp \
    generator.cpp \
    atomtype.cpp \
    moleculesystemcell.cpp \
    integrator/velocityverletintegrator.cpp \
    integrator/integrator.cpp \
    integrator/eulercromerintegrator.cpp \
    filemanager.cpp \
    modifier/modifier.cpp \
    modifier/berendsenthermostat.cpp \
    modifier/andersenthermostat.cpp \
    random.cpp \
    atom.cpp \
    force/lennardjonesforce.cpp \
    force/twoparticleforce.cpp \
    configurationparser.cpp

HEADERS += \
    moleculesystem.h \
    generator.h \
    atomtype.h \
    moleculesystemcell.h \
    integrator/velocityverletintegrator.h \
    integrator/integrator.h \
    integrator/eulercromerintegrator.h \
    filemanager.h \
    modifier/modifier.h \
    modifier/berendsenthermostat.h \
    modifier/andersenthermostat.h \
    random.h \
    atom.h \
    force/lennardjonesforce.h \
    force/twoparticleforce.h \
    configurationparser.h

# Building
myscript.target = myscript
myscript.commands = python $$PWD/../myscript.py $$PWD/../
QMAKE_EXTRA_TARGETS += myscript
PRE_TARGETDEPS += myscript

OTHER_FILES += ../testconfig.cfg \
    ../molecular-dynamics.pri
