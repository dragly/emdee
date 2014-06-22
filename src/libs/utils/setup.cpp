#include "setup.h"
#include "utils/logging.h"

/*!
 * \class Setup
 * \brief The Setup class is a helper class to load MPI and the glog library upon
 * starting an application.
 */

Setup::Setup(int argc, char* argv[])
#ifdef MD_USE_MPI
    :
    m_env(argc, argv)
  #endif
{
    setlocale(LC_ALL, "C");
#ifdef MD_USE_GLOG
    FLAGS_log_dir = ".";
    google::InitGoogleLogging(argv[0]);
#endif
}
