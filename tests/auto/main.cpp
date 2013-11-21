#include <unittest++/UnitTest++.h>
#include <fenv.h>

//#include <boost/mpi.hpp>
//namespace mpi = boost::mpi;

int main(int argc, char* argv[])
{
//    mpi::environment env(argc, argv);
//    mpi::communicator world;
    feenableexcept(FE_INVALID | FE_OVERFLOW);
    return UnitTest::RunAllTests();
}

