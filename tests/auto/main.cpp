#include <unittest++/UnitTest++.h>

#include <boost/mpi.hpp>
namespace mpi = boost::mpi;

int main(int argc, char* argv[])
{
    mpi::environment env(argc, argv);
    mpi::communicator world;
    return UnitTest::RunAllTests();
}

