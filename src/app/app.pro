TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

include(../../defaults.pri)

LIBS += -larmadillo -llapack -lblas -lconfig++

SOURCES += \
    main.cpp \
    configurationparser.cpp

HEADERS += \
    configurationparser.h

mpi {
    LIBS += -L../libs -lemdeempi
} else {
    LIBS += -L../libs -lemdee
}

OTHER_FILES += testconfig.cfg

copyconfig.commands = $(COPY_DIR) $$PWD/testconfig.cfg $$OUT_PWD
first.depends = $(first) copyconfig
export(first.depends)
export(copyconfig.commands)
QMAKE_EXTRA_TARGETS += first copyconfig
