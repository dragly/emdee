TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

include(../../defaults.pri)

mpi {
    LIBS += -L../../src/libs -lemdeempi
} else {
    LIBS += -L../../src/libs -lemdee
}

SOURCES += main.cpp

