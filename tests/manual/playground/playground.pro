#-------------------------------------------------
#
# Project created by QtCreator 2013-01-27T11:53:14
#
#-------------------------------------------------

include(../../../molecular-dynamics.pri)


QT       += testlib

QT       -= gui

TARGET = tst_playgroundtest
CONFIG   += console
CONFIG   -= app_bundle

TEMPLATE = app

INCLUDEPATH += $$SRC_DIR

SOURCES += tst_playgroundtest.cpp \
    ../../../src/moleculesystem.cpp \
    ../../../src/molecule.cpp \
    ../../../src/integrator.cpp \
    ../../../src/generator.cpp \
    ../../../src/atomtype.cpp \
    ../../../src/atom.cpp \
    ../../../src/moleculesystemcell.cpp \
    ../../../src/interatomicforce.cpp

DEFINES += SRCDIR=\\\"$$PWD/\\\"

HEADERS += \
    ../../../src/moleculesystem.h \
    ../../../src/molecule.h \
    ../../../src/integrator.h \
    ../../../src/generator.h \
    ../../../src/atomtype.h \
    ../../../src/atom.h \
    ../../../src/moleculesystemcell.h \
    ../../../src/interatomicforce.h
