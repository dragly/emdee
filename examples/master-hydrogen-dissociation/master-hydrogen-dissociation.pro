include(../../defaults.pri)
include(../../defaults-nonsrc.pri)

TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp

LIBS += -lyaml-cpp

OTHER_FILES += \
    configs/low_density_14K.yaml \
    configs/low_density_156K.yaml \
    configs/high_density_15600K.yaml \
    configs/high_density_156K.yaml \
    configs/high_density_14K.yaml \
    configs/low_density_15600K.yaml

