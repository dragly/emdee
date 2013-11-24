TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

include(../../defaults.pri)

LIBS += -lunittest++ -L../../src/libs -lemdee

SOURCES += main.cpp \
    forcestest.cpp \
    generatortest.cpp \
    systemtest.cpp \
    forcecelltest.cpp
