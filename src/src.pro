TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

TARGET = molecular-dynamics

SOURCES += main.cpp \
    moleculesystem.cpp \
    atom.cpp \
    generator.cpp \
    atomtype.cpp \
    molecule.cpp \
    integrator.cpp

HEADERS += \
    moleculesystem.h \
    atom.h \
    generator.h \
    atomtype.h \
    molecule.h \
    integrator.h

# Libraries
LIBS += -larmadillo -llapack -lblas

# Building
myscript.target = myscript
myscript.commands = python $$PWD/../myscript.py $$PWD/../
QMAKE_EXTRA_TARGETS += myscript
PRE_TARGETDEPS += myscript

COMMON_CXXFLAGS = -std=c++0x
QMAKE_CXXFLAGS += $$COMMON_CXXFLAGS
QMAKE_CXXFLAGS_RELEASE += $$COMMON_CXXFLAGS

QMAKE_CXXFLAGS_RELEASE += -O3
QMAKE_CXXFLAGS_RELEASE -= -O2
