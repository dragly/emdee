#include <src/moleculesystem.h>
#include <src/configurationparser.h>

#include <iostream>
#include <iomanip>
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
    cout << "Starting molecular-dynamics" << endl;
    cout << "I am processor number " << world.rank() << " of " << world.size() << endl;
    string configFileName = "testconfig.cfg";
    if(argc > 1) {
        configFileName = argv[1];
    }
    MoleculeSystem system;
    ConfigurationParser parser(&system);
    parser.runConfiguration(configFileName);
    cout << "Finished after " << setprecision(5) << timer.elapsed() << " seconds" << endl;
    return 0;
}

