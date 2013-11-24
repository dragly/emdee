TEMPLATE = app

include(../../defaults.pri)

LIBS += -larmadillo -llapack -lblas -lconfig++

SOURCES += \
    main.cpp \
    configurationparser.cpp

HEADERS += \
    configurationparser.h
