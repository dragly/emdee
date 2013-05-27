include(../../molecular-dynamics.pri)

TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

INCLUDEPATH = $$ROOT_DIR

SOURCES += $$system(find $$SRC_DIR -name \'*.cpp\')
SOURCES = $$replace(SOURCES, $$SRC_DIR/main.cpp, )

SOURCES += main.cpp

DEFINES += ARMA_NO_DEBUG
