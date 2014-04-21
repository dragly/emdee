CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

include(../../defaults.pri)

LIBS += -lunittest++
mpi {
    LIBS += -L../../src/libs -lemdeempi
} else {
    LIBS += -L../../src/libs -lemdee
}

SOURCES += main.cpp \
    forcestest.cpp \
    generatortest.cpp \
    systemtest.cpp \
    forcecelltest.cpp \
    threeparticleforcesystem.cpp \
    fanntest.cpp \
    twoparticleforcecount.cpp
