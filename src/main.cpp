#include <src/moleculesystem.h>
#include <src/configurationparser.h>

#include <iostream>
#include <libconfig.h++>

//#include <openmpi/mpi.h>

using namespace std;
using namespace libconfig;

int main(int argc, char** argv)
{
//    MPI_Init(&argc, &argv);
    string configFileName = "testconfig.cfg";
    if(argc > 1) {
        configFileName = argv[1];
    }
    MoleculeSystem system;
    ConfigurationParser parser(&system);
    parser.runConfiguration(configFileName);
//    MPI_Finalize();
    return 0;
}

