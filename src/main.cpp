#include <src/moleculesystem.h>
#include <src/configurationparser.h>

#include <iostream>
#include <libconfig.h++>

#include <boost/mpi.hpp>

//#include <openmpi/mpi.h>

using namespace std;
using namespace libconfig;
namespace mpi = boost::mpi;

int main(int argc, char** argv)
{
    mpi::environment env(argc, argv);
    mpi::communicator world;
    mpi::timer timer;
    timer.restart();
    cout << "I am processor number " << world.rank() << endl;
    string configFileName = "testconfig.cfg";
    if(argc > 1) {
        configFileName = argv[1];
    }
    MoleculeSystem system;
    ConfigurationParser parser(&system);
    parser.runConfiguration(configFileName);
    cout << "Finished after " << timer.elapsed() << " seconds" << endl;
    return 0;
}

