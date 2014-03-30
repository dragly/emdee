#include <unittest++/UnitTest++.h>
#include <unittest++/Test.h>
#include <unittest++/TestReporterStdout.h>
#include <unittest++/TestRunner.h>

#include <fenv.h>
#include <utils/logging.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#ifdef MD_USE_MPI
#include <boost/mpi.hpp>
namespace mpi = boost::mpi;
#endif

int main(int argc, char* argv[])
{
#ifdef MD_USE_GLOG
    FLAGS_log_dir = ".";
    google::InitGoogleLogging(argv[0]);
#endif
#ifdef MD_USE_MPI
    mpi::environment env(argc, argv);
    mpi::communicator world;
    (void)env;
    (void)world;

#else
    (void)argc;
#endif
    feenableexcept(FE_INVALID | FE_OVERFLOW);

//    int result = 0;
//    UnitTest::TestReporterStdout reporter;
//    UnitTest::TestRunner runner(reporter);
//    result = runner.RunTestsIf(UnitTest::Test::GetTestList(), "ForceCellTest", UnitTest::True(), 0);
//    return result;

    return UnitTest::RunAllTests();
}

