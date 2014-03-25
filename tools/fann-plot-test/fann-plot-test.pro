TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt
CONFIG += c++11

SOURCES += main.cpp

LIBS += -lfann

QMAKE_CXX = ccache mpicxx

QMAKE_CXXFLAGS += $$system(mpicxx --showme:compile)
LIBS += $$system(mpicxx --showme:link)
