#ifndef PROCESSOR_H
#define PROCESSOR_H

#include <src/range.h>

#include <boost/mpi.hpp>
#include <vector>
#include <armadillo>
using namespace std;
using namespace arma;
namespace mpi = boost::mpi;

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

    const vector<MoleculeSystemCell*> cells() const;
    bool shouldSendFirst(const irowvec &direction);
protected:

    void receiveAtomsFromNeighbor(const ProcessorNeighbor &neighbor);
    void sendAtomsToNeighbor(const ProcessorNeighbor &neighbor);

    MoleculeSystem* m_moleculeSystem;

    mpi::communicator world;
    mpi::timer communicationTimer;
    mpi::timer pureCommunicationTimer;

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

    vector<MoleculeSystemCell*> m_cells;

    double m_totalCommunicationTime;
    double m_pureCommunicationTime;

    vector<ProcessorNeighbor> sendNeighbors;
    vector<ProcessorNeighbor> receiveNeighbors;

    vector<irowvec> directions;
};

inline const vector<MoleculeSystemCell*> Processor::cells() const {
    return m_cells;
}

#endif // PROCESSOR_H
