#include <moleculesystem.h>
#include <configurationparser.h>
#include <utils/logging.h>

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
#ifdef MD_USE_GLOG
    FLAGS_log_dir = ".";
    google::InitGoogleLogging(argv[0]);
#endif
    LOG(INFO) << "Starting molecular-dynamics";
#ifdef MD_USE_MPI
    mpi::environment env(argc, argv);
    mpi::communicator world;
    mpi::timer timer;
#endif
    feenableexcept(FE_INVALID | FE_OVERFLOW);
#ifdef MD_USE_MPI
    timer.restart();
    LOG(INFO) << "I am processor number " << world.rank() << " of " << world.size();
#endif
    string configFileName = "testconfig.cfg";
    if(argc > 1 && argv[1][0] != '-') {
        configFileName = argv[1];
    }
    MoleculeSystem system;
    system.setPeriodicity(true, true, true);
    ConfigurationParser parser(&system);
    parser.runConfiguration(configFileName);
#ifdef MD_USE_MPI
    LOG(INFO) << "Finished after " << setprecision(5) << timer.elapsed() << " seconds";
#endif


    return 0;
}

