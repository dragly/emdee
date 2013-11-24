#include <moleculesystem.h>
#include <configurationparser.h>

#include <iostream>
#include <iomanip>
#include <libconfig.h++>
#include <fenv.h>

//#include <boost/mpi.hpp>

//#include <openmpi/mpi.h>

using namespace std;
using namespace libconfig;
//namespace mpi = boost::mpi;

int main(int argc, char** argv)
{
    cout << "Starting molecular-dynamics" << endl;
#ifdef MD_USE_MPI
    mpi::environment env(argc, argv);
    mpi::communicator world;
    mpi::timer timer;
#endif
    feenableexcept(FE_INVALID | FE_OVERFLOW);
#ifdef MD_USE_MPI
    timer.restart();
    cout << "I am processor number " << world.rank() << " of " << world.size() << endl;
#endif
    string configFileName = "testconfig.cfg";
    if(argc > 1) {
        configFileName = argv[1];
    }
    MoleculeSystem system;
    ConfigurationParser parser(&system);
    parser.runConfiguration(configFileName);
#ifdef MD_USE_MPI
    cout << "Finished after " << setprecision(5) << timer.elapsed() << " seconds" << endl;
#endif
    return 0;
}

