include(../../molecular-dynamics.pri)

TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

INCLUDEPATH = $$ROOT_DIR

SOURCES += main.cpp \
    ../../src/interatomicforce.cpp \
    ../../src/generator.cpp \
    ../../src/filemanager.cpp \
    ../../src/atomtype.cpp \
    ../../src/atom.cpp \
    ../../src/random.cpp \
    ../../src/moleculesystemcell.cpp \
    ../../src/moleculesystem.cpp \
    ../../src/molecule.cpp \
    ../../src/integrator/velocityverletintegrator.cpp \
    ../../src/integrator/integrator.cpp \
    ../../src/integrator/eulercromerintegrator.cpp \
    ../../src/modifier/modifier.cpp \
    ../../src/modifier/berendsenthermostat.cpp \
    ../../src/modifier/andersenthermostat.cpp

HEADERS += \
    ../../src/interatomicforce.h \
    ../../src/generator.h \
    ../../src/filemanager.h \
    ../../src/atomtype.h \
    ../../src/atom.h \
    ../../src/random.h \
    ../../src/moleculesystemcell.h \
    ../../src/moleculesystem.h \
    ../../src/molecule.h \
    ../../src/integrator/velocityverletintegrator.h \
    ../../src/integrator/integrator.h \
    ../../src/integrator/eulercromerintegrator.h \
    ../../src/modifier/modifier.h \
    ../../src/modifier/berendsenthermostat.h \
    ../../src/modifier/andersenthermostat.h
