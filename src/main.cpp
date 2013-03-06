#include <src/moleculesystem.h>
#include <src/configurationparser.h>

#include <iostream>
#include <libconfig.h++>

#include <boost/mpi.hpp>

using namespace std;
using namespace libconfig;
namespace mpi = boost::mpi;

int main(int argc, char** argv)
{
//    MPI_Init(&argc, &argv);
    mpi::timer totalTimer;
    totalTimer.restart();
    string configFileName = "testconfig.cfg";
    if(argc > 1) {
        configFileName = argv[1];
    }
    MoleculeSystem system;
    ConfigurationParser parser(&system);
    parser.runConfiguration(configFileName);
//    MPI_Finalize();
    cout << "Total run time " << totalTimer.elapsed() << endl;
    return 0;
}

