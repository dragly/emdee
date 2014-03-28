CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

include(../../defaults.pri)

LIBS += -lunittest++
!nompi {
    LIBS += -L../../src/libs -lemdeempi
} else {
    LIBS += -L../../src/libs -lemdee
}

SOURCES += main.cpp \
    forcestest.cpp \
    generatortest.cpp \
    systemtest.cpp \
    forcecelltest.cpp
