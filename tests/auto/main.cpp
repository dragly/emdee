#include <unittest++/UnitTest++.h>
#include <fenv.h>
#include <glog/logging.h>

#ifdef MD_USE_MPI
#include <boost/mpi.hpp>
namespace mpi = boost::mpi;
#endif

int main(int argc, char* argv[])
{
    FLAGS_log_dir = ".";
    google::InitGoogleLogging(argv[0]);
#ifdef MD_USE_MPI
    mpi::environment env(argc, argv);
    mpi::communicator world;
#else
    (void)argc;
#endif
    feenableexcept(FE_INVALID | FE_OVERFLOW);
    return UnitTest::RunAllTests();
}

