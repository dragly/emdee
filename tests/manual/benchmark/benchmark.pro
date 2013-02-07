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
    ../../../src/molecule.cpp \
    ../../../src/interatomicforce.cpp \
    ../../../src/integrator.cpp \
    ../../../src/generator.cpp \
    ../../../src/atomtype.cpp \
    ../../../src/atom.cpp

HEADERS += \
    ../../../src/moleculesystemcell.h \
    ../../../src/moleculesystem.h \
    ../../../src/molecule.h \
    ../../../src/interatomicforce.h \
    ../../../src/integrator.h \
    ../../../src/generator.h \
    ../../../src/atomtype.h \
    ../../../src/atom.h
