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
    configurationparser.cpp \
    math/vector3.cpp \
    processor.cpp \
    range.cpp \
    modifier/constantforce.cpp

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
    configurationparser.h \
    math/vector3.h \
    processor.h \
    range.h \
    modifier/constantforce.h \
    progressreporter.h

#CURRENT_BRANCH = $$system(git rev-parse --abbrev-ref HEAD)

#OBJECTS_DIR = $$CURRENT_BRANCH

#system(mkdir -p $$CURRENT_BRANCH)

# Building
myscript.target = myscript
myscript.commands = python $$PWD/../myscript.py $$PWD/../ $$CURRENT_BRANCH
QMAKE_EXTRA_TARGETS += myscript
PRE_TARGETDEPS += myscript

OTHER_FILES += ../testconfig.cfg \
    ../molecular-dynamics.pri
