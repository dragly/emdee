#-------------------------------------------------
#
# Project created by QtCreator 2013-02-11T13:08:27
#
#-------------------------------------------------

include(../../tests.pri)

QT       += testlib

QT       -= gui

TARGET = tst_generatortest
CONFIG   += console
CONFIG   -= app_bundle

TEMPLATE = app


SOURCES += tst_generatortest.cpp \
    ../../../src/generator.cpp \
    ../../../src/atomtype.cpp \
    ../../../src/atom.cpp \
    ../../../src/vector3d.cpp
DEFINES += SRCDIR=\\\"$$PWD/\\\"

HEADERS += \
    ../../../src/generator.h \
    ../../../src/atomtype.h \
    ../../../src/atom.h \
    ../../../src/vector3d.h
