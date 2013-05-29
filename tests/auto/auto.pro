TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

include(../../molecular-dynamics.pri)

LIBS += -lunittest++

SOURCES += main.cpp \
    forcestest.cpp \
    generatortest.cpp \
    systemtest.cpp \
    forcecelltest.cpp
message(find $$SRC_DIR -name '*.cpp')
SOURCES += $$system(find $$SRC_DIR -name \'*.cpp\')
SOURCES = $$replace(SOURCES, $$SRC_DIR/main.cpp, )

message($$SOURCES)
