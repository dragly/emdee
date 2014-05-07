ROOT_DIR = $$PWD
SRC_DIR = $$PWD/src
INCLUDEPATH += $$ROOT_DIR/src/libs

message($$INCLUDEPATH)

COMMON_CXXFLAGS = -std=c++11
QMAKE_CXXFLAGS += $$COMMON_CXXFLAGS
QMAKE_CXXFLAGS_RELEASE += $$COMMON_CXXFLAGS
QMAKE_CXXFLAGS_DEBUG += $$COMMON_CXXFLAGS

CONFIG(release, debug|release) {
    DEFINES += NDEBUG # used by Google logging library
    DEFINES += ARMA_NO_DEBUG
    QMAKE_CXXFLAGS_RELEASE -= -O2
    QMAKE_CXXFLAGS_RELEASE += -O3
}

android {
    INCLUDEPATH += $$ROOT_DIR/armadillo
}

LIBS += -lglog -ldoublefann

!noglog {
    DEFINES += MD_USE_GLOG
}

mpi {
    DEFINES += MD_USE_MPI
    DEFINES+=USE_MPI

    QMAKE_CXX = mpicxx
    QMAKE_CXX_RELEASE = $$QMAKE_CXX
    QMAKE_CXX_DEBUG = $$QMAKE_CXX
    QMAKE_LINK = $$QMAKE_CXX

    QMAKE_CFLAGS += $$system(mpicc --showme:compile)
    QMAKE_LFLAGS += $$system(mpicxx --showme:link)
    QMAKE_CXXFLAGS += $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK
    QMAKE_CXXFLAGS_RELEASE += $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK

    LIBS += -lboost_system -lboost_filesystem -lboost_mpi -lboost_serialization
}

QMAKE_CXXFLAGS += -g
QMAKE_CXXFLAGS_RELEASE += -g
# -ftree-parallelize-loops=4 -ftree-loop-optimize -floop-parallelize-all

dev {
    DEFINES += DEVELOPMENT_TESTS
}

#!noccache {
#    QMAKE_CXX = ccache $$QMAKE_CXX
#}
