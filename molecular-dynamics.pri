ROOT_DIR = $$PWD
SRC_DIR = $$PWD/src

# Libraries
LIBS += -larmadillo -llapack -lblas -lconfig++ -lhdf5_cpp -lhdf5

COMMON_CXXFLAGS = -std=c++0x
QMAKE_CXXFLAGS += $$COMMON_CXXFLAGS
QMAKE_CXXFLAGS_RELEASE += $$COMMON_CXXFLAGS

QMAKE_CXXFLAGS_RELEASE += -O3
QMAKE_CXXFLAGS_RELEASE -= -O2

release {
    DEFINES += ARMA_NO_DEBUG
}
