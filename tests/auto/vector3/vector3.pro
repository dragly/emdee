#-------------------------------------------------
#
# Project created by QtCreator 2013-03-01T14:17:52
#
#-------------------------------------------------

include(../../tests.pri)

QT       += testlib

QT       -= gui

TARGET = tst_vector3test
CONFIG   += console
CONFIG   -= app_bundle

TEMPLATE = app

INCLUDEPATH += $$ROOT_DIR


SOURCES += tst_vector3test.cpp \
    ../../../src/math/vector3.cpp
DEFINES += SRCDIR=\\\"$$PWD/\\\"

HEADERS += \
    ../../../src/math/vector3.h
