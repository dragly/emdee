#ifndef PROCESSOR_H
#define PROCESSOR_H

#include <range.h>
#include <vector>
#include <armadillo>
using namespace std;
using namespace arma;

class MoleculeSystem;
class MoleculeSystemCell;
class Atom;

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
//    void receiveForcesFromNeighbor(const ProcessorNeighbor &neighbor);
//    void sendForcesToNeighbor(const ProcessorNeighbor &neighbor);
    //    void communicateForces();
    void communicateForces();
//    bool checkDirection(const irowvec &direction);
    void clearForcesInNeighborCells();
protected:

    MoleculeSystem* m_moleculeSystem;

#ifdef USE_MPI
    mpi::communicator world;
    mpi::timer communicationTimer;
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

    vector<irowvec> directions;
    vector<irowvec> forceDirections;
};

#endif // PROCESSOR_H
