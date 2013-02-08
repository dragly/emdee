#-------------------------------------------------
#
# Project created by QtCreator 2013-02-07T13:26:13
#
#-------------------------------------------------

include(../../tests.pri)

QT       += testlib

QT       -= gui

TARGET = tst_systemtest
CONFIG   += console
CONFIG   -= app_bundle

TEMPLATE = app


SOURCES += tst_systemtest.cpp \
    ../../../src/moleculesystemcell.cpp \
    ../../../src/moleculesystem.cpp \
    ../../../src/molecule.cpp \
    ../../../src/interatomicforce.cpp \
    ../../../src/generator.cpp \
    ../../../src/atomtype.cpp \
    ../../../src/atom.cpp \
    ../../../src/integrator/velocityverletintegrator.cpp \
    ../../../src/integrator/integrator.cpp \
    ../../../src/integrator/eulercromerintegrator.cpp

DEFINES += SRCDIR=\\\"$$PWD/\\\"

HEADERS += \
    ../../../src/moleculesystemcell.h \
    ../../../src/moleculesystem.h \
    ../../../src/molecule.h \
    ../../../src/interatomicforce.h \
    ../../../src/generator.h \
    ../../../src/atomtype.h \
    ../../../src/atom.h \
    ../../../src/integrator/velocityverletintegrator.h \
    ../../../src/integrator/integrator.h \
    ../../../src/integrator/eulercromerintegrator.h
