#-------------------------------------------------
#
# Project created by QtCreator 2013-02-06T13:46:28
#
#-------------------------------------------------

include(../../tests.pri)

QT       += testlib

QT       -= gui

TARGET = tst_forcestest
CONFIG   += console
CONFIG   -= app_bundle

TEMPLATE = app

INCLUDEPATH += ../../../src

SOURCES += tst_forcestest.cpp \
    ../../../src/interatomicforce.cpp \
    ../../../src/atomtype.cpp \
    ../../../src/atom.cpp
DEFINES += SRCDIR=\\\"$$PWD/\\\"

HEADERS += \
    ../../../src/interatomicforce.h \
    ../../../src/atomtype.h \
    ../../../src/atom.h
