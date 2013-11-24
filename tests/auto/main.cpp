#include <unittest++/UnitTest++.h>
#include <fenv.h>

#ifdef MD_USE_MPI
#include <boost/mpi.hpp>
namespace mpi = boost::mpi;
#endif

int main(int argc, char* argv[])
{
#ifdef MD_USE_MPI
    mpi::environment env(argc, argv);
    mpi::communicator world;
#endif
    feenableexcept(FE_INVALID | FE_OVERFLOW);
    return UnitTest::RunAllTests();
}

