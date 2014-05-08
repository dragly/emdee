#ifndef SETUP_H
#define SETUP_H

#ifdef MD_USE_MPI
#include <boost/mpi.hpp>
#endif

class Setup
{
public:
    Setup(int argc, char *argv[]);
private:
#ifdef MD_USE_MPI
    boost::mpi::environment m_env;
    boost::mpi::communicator m_world;
#endif
};

#endif // SETUP_H
