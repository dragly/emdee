#-------------------------------------------------
#
# Project created by QtCreator 2013-02-07T09:31:13
#
#-------------------------------------------------

include(../../tests.pri)

QT       += testlib

QT       -= gui

TARGET = tst_benchmarktest
CONFIG   += console
CONFIG   -= app_bundle

TEMPLATE = app

DEFINES += SRCDIR=\\\"$$PWD/\\\"


SOURCES += tst_benchmarktest.cpp \
    ../../../src/moleculesystemcell.cpp \
    ../../../src/moleculesystem.cpp \
    ../../../src/force/twoparticleforce.cpp \
    ../../../src/generator.cpp \
    ../../../src/atomtype.cpp \
    ../../../src/atom.cpp \
    ../../../src/integrator/velocityverletintegrator.cpp \
    ../../../src/integrator/integrator.cpp \
    ../../../src/integrator/eulercromerintegrator.cpp \
    ../../../src/filemanager.cpp \
    ../../../src/force/lennardjonesforce.cpp \
    ../../../src/math/vector3.cpp \
    ../../../src/processor.cpp \
    ../../../src/range.cpp

HEADERS += \
    ../../../src/moleculesystemcell.h \
    ../../../src/moleculesystem.h \
    ../../../src/force/twoparticleforce.h \
    ../../../src/generator.h \
    ../../../src/atomtype.h \
    ../../../src/atom.h \
    ../../../src/integrator/velocityverletintegrator.h \
    ../../../src/integrator/integrator.h \
    ../../../src/integrator/eulercromerintegrator.h \
    ../../../src/filemanager.h \
    ../../../src/force/lennardjonesforce.h \
    ../../../src/math/vector3.h \
    ../../../src/processor.h \
    ../../../src/range.h
