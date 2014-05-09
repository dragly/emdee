TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

include(../../defaults.pri)
include(../../defaults-nonsrc.pri)

LIBS += -lconfig++

SOURCES += \
    main.cpp \
    configurationparser.cpp

HEADERS += \
    configurationparser.h

OTHER_FILES += testconfig.cfg

copyconfig.commands = $(COPY_DIR) $$PWD/testconfig.cfg $$OUT_PWD
first.depends = $(first) copyconfig
export(first.depends)
export(copyconfig.commands)
QMAKE_EXTRA_TARGETS += first copyconfig
