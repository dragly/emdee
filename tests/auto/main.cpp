#include <unittest++/UnitTest++.h>
#include <unittest++/Test.h>
#include <unittest++/TestReporterStdout.h>
#include <unittest++/TestRunner.h>

#include <fenv.h>
#include <utils/logging.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <locale>
#include <vector>
#include <string>

#ifdef MD_USE_MPI
#include <boost/mpi.hpp>
namespace mpi = boost::mpi;
#endif

int main(int argc, char* argv[])
{
    setlocale(LC_ALL, "C");
#ifdef MD_USE_GLOG
    FLAGS_log_dir = ".";
    google::InitGoogleLogging(argv[0]);
#endif
#ifdef MD_USE_MPI
    mpi::environment env(argc, argv);
    mpi::communicator world;
    (void)env;
    (void)world;
#endif
    feenableexcept(FE_INVALID | FE_OVERFLOW);


    std::vector<std::string> tests;
    if(argc > 1 && !strcmp(argv[1], "dev")) {
        tests.push_back("FannForceSystem");
    } else if(argc > 1 && !strcmp(argv[1], "dev2")) {
        tests.push_back("Development2");
    } else {
        tests.push_back("FannDerivative");
        tests.push_back("ForceCell");
        tests.push_back("Forces");
        tests.push_back("System");
        tests.push_back("Generator");
        tests.push_back("ThreeParticleForceSystem");
        tests.push_back("TwoParticleForceCount");
    }
    int result = 0;
    for(const std::string& testName : tests) {
        std::cout << "Running " << testName << std::endl;
        UnitTest::TestReporterStdout reporter;
        UnitTest::TestRunner runner(reporter);
        result += runner.RunTestsIf(UnitTest::Test::GetTestList(), testName.c_str(), UnitTest::True(), 0);
    }
    if(result != 0) {
        std::cerr << "FAILURE: Some tests failed! See above log for details..." << std::endl;
    }
    return result;
    //    return UnitTest::RunAllTests();
}

