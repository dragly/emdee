#ifndef PROCESSOR_H
#define PROCESSOR_H

#include <range.h>
#ifdef USE_MPI
#include <boost/mpi.hpp>
namespace mpi = boost::mpi;
#endif
#include <vector>
#include <armadillo>
using namespace std;
using namespace arma;

class MoleculeSystem;
class MoleculeSystemCell;
class Atom;

class ProcessorNeighbor {
public:
    ProcessorNeighbor();

    vector<MoleculeSystemCell*> cells;
//    vector<MoleculeSystemCell*> cellsToReceive;
    irowvec direction;
    irowvec coordinates;
    int rank;
};

class Processor
{
public:
    Processor(MoleculeSystem* moleculeSystem);

    void setupProcessors();

    int rank();
    void communicateAtoms();

    const vector<MoleculeSystemCell *> &localCells() const;
    const vector<MoleculeSystemCell *> &localAndGhostCells() const;
    bool shouldSendFirst(const irowvec &direction);
    int nAtoms();
    int nProcessors();
    void clearForcesInNeighborCells();
protected:

    void receiveAtomsFromNeighbor(const ProcessorNeighbor &neighbor);
    void sendAtomsToNeighbor(const ProcessorNeighbor &neighbor);

    MoleculeSystem* m_moleculeSystem;

#ifdef USE_MPI
    mpi::communicator world;
    mpi::timer m_communicationTimer;
    mpi::timer pureCommunicationTimer;
#endif

    int m_nProcessors;

    int m_nProcessorsX;
    int m_nProcessorsY;
    int m_nProcessorsZ;

    int m_coordinateX;
    int m_coordinateY;
    int m_coordinateZ;

    Range m_cellRangeX;
    Range m_cellRangeY;
    Range m_cellRangeZ;

    vector<MoleculeSystemCell*> m_localCells;
    vector<MoleculeSystemCell*> m_localAndGhostCells;

    double m_totalCommunicationTime;
    double m_pureCommunicationTime;

    vector<ProcessorNeighbor> m_sendNeighbors;
    vector<ProcessorNeighbor> m_receiveNeighbors;

    vector<ProcessorNeighbor> forceSendNeighbors;
    vector<ProcessorNeighbor> forceReceiveNeighbors;

    vector<irowvec> m_directions;
    vector<irowvec> forceDirections;
};

inline const vector<MoleculeSystemCell *> &Processor::localCells() const {
    return m_localCells;
}

inline const vector<MoleculeSystemCell *> &Processor::localAndGhostCells() const
{
    return m_localAndGhostCells;
}

#endif // PROCESSOR_H
