TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

include(../../molecular-dynamics.pri)

LIBS += -lunittest++

SOURCES += main.cpp \
    forcestest.cpp \
    generatortest.cpp \
    systemtest.cpp

SOURCES += $$SRC_DIR/*.cpp \
            $$SRC_DIR/force/*.cpp \
            $$SRC_DIR/integrator/*.cpp \
            $$SRC_DIR/math/*.cpp \
            $$SRC_DIR/modifier/*.cpp

SOURCES = $$system(ls $$SOURCES)
SOURCES = $$replace(SOURCES, $$SRC_DIR/main.cpp, )
