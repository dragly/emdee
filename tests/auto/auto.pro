TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

include(../../molecular-dynamics.pri)

LIBS += -lunittest++

SOURCES += main.cpp \
    forcestest.cpp \
    ../../src/force/twoparticleforce.cpp \
    ../../src/force/singleparticleforce.cpp \
    ../../src/force/lennardjonesforce.cpp \
    ../../src/force/constantforce.cpp \
    ../../src/integrator/velocityverletintegrator.cpp \
    ../../src/integrator/integrator.cpp \
    ../../src/integrator/eulercromerintegrator.cpp

HEADERS += \
    ../../src/force/twoparticleforce.h \
    ../../src/force/singleparticleforce.h \
    ../../src/force/lennardjonesforce.h \
    ../../src/force/constantforce.h \
    ../../src/integrator/velocityverletintegrator.h \
    ../../src/integrator/integrator.h \
    ../../src/integrator/eulercromerintegrator.h

