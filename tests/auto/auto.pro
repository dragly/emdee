CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

include(../../defaults.pri)
include(../../defaults-nonsrc.pri)

LIBS += -lunittest++

SOURCES += main.cpp \
    forcestest.cpp \
    generatortest.cpp \
    systemtest.cpp \
    forcecelltest.cpp \
    threeparticleforcesystem.cpp \
    fanntest.cpp \
    twoparticleforcecount.cpp \
    fannderivative.cpp \
    development2.cpp
