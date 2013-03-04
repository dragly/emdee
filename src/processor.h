#ifndef PROCESSOR_H
#define PROCESSOR_H

#include <src/range.h>

#include <boost/mpi.hpp>
#include <vector>
using namespace std;
namespace mpi = boost::mpi;

class MoleculeSystem;
class MoleculeSystemCell;

class Processor
{
public:
    Processor(MoleculeSystem* moleculeSystem);

    void setupProcessors();

    int rank();
protected:
    MoleculeSystem* m_moleculeSystem;

    mpi::communicator world;

    int m_nProcessors;

    int m_nProcessorsX;
    int m_nProcessorsY;
    int m_nProcessorsZ;

    int m_coordinateX;
    int m_coordinateY;
    int m_coordinateZ;

    Range m_rangeX;
    Range m_rangeY;
    Range m_rangeZ;

    vector<MoleculeSystemCell*> m_cells;
};

#endif // PROCESSOR_H
