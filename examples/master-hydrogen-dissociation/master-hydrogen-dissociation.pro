include(../../defaults.pri)
include(../../defaults-nonsrc.pri)

TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp

LIBS += -lyaml-cpp

OTHER_FILES += \
    configs/low_density_14K.yaml

