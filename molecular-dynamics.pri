ROOT_DIR = $$PWD
SRC_DIR = $$PWD/src
QMAKE_CXX = ccache g++
# Libraries
LIBS += -larmadillo -llapack -lblas -lconfig++ -lhdf5_cpp -lhdf5 -lboost_filesystem -lboost_system

COMMON_CXXFLAGS = -std=c++0x
QMAKE_CXXFLAGS += $$COMMON_CXXFLAGS
QMAKE_CXXFLAGS_RELEASE += $$COMMON_CXXFLAGS
QMAKE_CXXFLAGS_DEBUG += $$COMMON_CXXFLAGS

release {
    DEFINES += ARMA_NO_DEBUG
    QMAKE_CXXFLAGS_RELEASE -= -O2
    QMAKE_CXXFLAGS_RELEASE += -O3
}
CONFIG += mpi
mpi {
    QMAKE_CXX = ccache mpicxx
    QMAKE_CXX_RELEASE = $$QMAKE_CXX
    QMAKE_CXX_DEBUG = $$QMAKE_CXX
    QMAKE_LINK = $$QMAKE_CXX

    QMAKE_CFLAGS += $$system(mpicc --showme:compile)
    QMAKE_LFLAGS += $$system(mpicxx --showme:link)
    QMAKE_CXXFLAGS += $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK
    QMAKE_CXXFLAGS_RELEASE += $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK
}

# -ftree-parallelize-loops=4 -ftree-loop-optimize -floop-parallelize-all
