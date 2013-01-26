TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

myscript.target = myscript
myscript.commands = python $$PWD/myscript.py $$PWD
QMAKE_EXTRA_TARGETS += myscript
PRE_TARGETDEPS += myscript

SOURCES += main.cpp \
    moleculesystem.cpp \
    atom.cpp

HEADERS += \
    moleculesystem.h \
    atom.h

