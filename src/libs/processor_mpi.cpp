#include "processor.h"
#include <moleculesystem.h>
#include <moleculesystemcell.h>
#include <range.h>
#include <atom.h>
#include <math/vector3.h>
#include <utils/logging.h>

#include <iostream>
#include <iomanip>
#include <math.h>

#include <boost/mpi.hpp>
#include <boost/mpi/collectives.hpp>
namespace mpi = boost::mpi;

using namespace std;

Processor::Processor(MoleculeSystem *moleculeSystem) :
    m_moleculeSystem(moleculeSystem),
    m_nProcessors(0),
    m_nProcessorsX(0),
    m_nProcessorsY(0),
    m_nProcessorsZ(0),
    m_totalCommunicationTime(0),
    m_pureCommunicationTime(0)
{
    LOG(INFO) << "Processor loaded. I am rank " << world.rank() << " of " << world.size() << " processors.";
    irowvec left = {-1, 0, 0};
    m_directions.push_back(left);
    irowvec right = {1, 0, 0};
    m_directions.push_back(right);
    irowvec front = {0, 1, 0};
    m_directions.push_back(front);
    irowvec back = {0, -1, 0};
    m_directions.push_back(back);
    irowvec down = {0, 0, -1};
    m_directions.push_back(down);
    irowvec up = {0, 0, 1};
    m_directions.push_back(up);

    for(int i = -1; i < 2; i++) {
        for(int j = -1; j < 2; j++) {
            for(int k = -1; k < 2; k++) {
                irowvec direction = {i, j, k};
                if(
                        !( // if not one of ..
                           ((direction(0) >= 0 && direction(1) >= 0) && !(direction(0) == 0 && direction(1) == 0 && direction(2) == -1)) // 2x2 in upper right (except right down)
                           || (direction(0) == 1 && direction(1) == -1)) // 1x1 lower right
                        ) // then continue ...
                    //                neighbor->hasAlreadyCalculatedForcesBetweenSelfAndNeighbors()
                    //            )
                {
                    continue;
                }
                forceDirections.push_back(direction);
            }
        }
    }
}

void Processor::setupProcessors()
{
    m_sendNeighbors.clear();
    m_receiveNeighbors.clear();
    m_localCells.clear();
    m_localAndGhostCells.clear();

    m_nProcessors = world.size();
    //    cout << "Setting up " << m_nProcessors << " processors" << endl;
    m_nProcessorsZ = int( pow(m_nProcessors, 1./3.) + 1e-9);
    m_nProcessorsY = int( sqrt((double)m_nProcessors / (double)m_nProcessorsZ) + 1e-9);
    m_nProcessorsX = int( (double)m_nProcessors / ((double)m_nProcessorsZ * (double) m_nProcessorsY) + 1e-9);
    //    cout << "Dividing processor space into " << m_nProcessorsX << " * " << m_nProcessorsY << " * " << m_nProcessorsZ << " processors.";
    //    cout << "This gives " << m_moleculeSystem->nCells()(0) / (double)m_nProcessorsX << " * " << m_moleculeSystem->nCells()(1) / (double)m_nProcessorsY << " * " << m_moleculeSystem->nCells()(2) / (double)m_nProcessorsZ << " cells per processor" << endl;

    // Find own processor coordinates
    m_coordinateX = world.rank() / (m_nProcessorsZ * m_nProcessorsY);
    m_coordinateY = (world.rank() / m_nProcessorsZ) % m_nProcessorsY;
    m_coordinateZ = world.rank() % m_nProcessorsZ;

//        cout << "My coordinate is " << m_coordinateX << " " << m_coordinateY << " " << m_coordinateZ << endl;

    // Find my cells
    m_cellRangeX = Range(m_coordinateX, m_nProcessorsX, m_moleculeSystem->nCells()(0));
    m_cellRangeY = Range(m_coordinateY, m_nProcessorsY, m_moleculeSystem->nCells()(1));
    m_cellRangeZ = Range(m_coordinateZ, m_nProcessorsZ, m_moleculeSystem->nCells()(2));

//        cout << "X Range: " << m_cellRangeX << endl;
//        cout << "Y Range: " << m_cellRangeY << endl;
//        cout << "Z Range: " << m_cellRangeZ << endl;

    // Clear local cell state
    for(MoleculeSystemCell* cell : m_moleculeSystem->globalCells()) {
        cell->setLocal(false);
    }

    // Find all local cells
    for(int i = m_cellRangeX.firstElement(); i <= m_cellRangeX.lastElement(); i++) {
        for(int j = m_cellRangeY.firstElement(); j <= m_cellRangeY.lastElement(); j++) {
            for(int k = m_cellRangeZ.firstElement(); k <= m_cellRangeZ.lastElement(); k++) {
                //                cout << "I have cells " << i << " " << j << " " << k << endl;
                MoleculeSystemCell* cell = m_moleculeSystem->cell(i,j,k);
//                if((m_nProcessorsX > 1 && m_cellRangeX.firstElement() == i) || (m_nProcessorsY > 1 && m_cellRangeY.firstElement() == j) || (m_nProcessorsZ > 1 && m_cellRangeZ.firstElement() == k) ||
//                        (m_nProcessorsX > 1 && m_cellRangeX.lastElement() == i) || (m_nProcessorsY > 1 && m_cellRangeY.lastElement() == j) || (m_nProcessorsZ > 1 && m_cellRangeZ.lastElement() == k)) {
//                    cell->setOnProcessorEdge(true);
//                }
                cell->setLocal(true);
                m_localCells.push_back(cell);
                m_localAndGhostCells.push_back(cell);
            }
        }
    }

    // Find neighbors
    for(const irowvec& direction : m_directions) {

        ProcessorNeighbor sendNeighbor;
        ProcessorNeighbor receiveNeighbor;

        sendNeighbor.direction = direction;
        receiveNeighbor.direction = -direction;

        sendNeighbor.coordinates(0) = (m_coordinateX + sendNeighbor.direction(0) + m_nProcessorsX) % m_nProcessorsX;
        sendNeighbor.coordinates(1) = (m_coordinateY + sendNeighbor.direction(1) + m_nProcessorsY) % m_nProcessorsY;
        sendNeighbor.coordinates(2) = (m_coordinateZ + sendNeighbor.direction(2) + m_nProcessorsZ) % m_nProcessorsZ;

        receiveNeighbor.coordinates(0) = (m_coordinateX + receiveNeighbor.direction(0) + m_nProcessorsX) % m_nProcessorsX;
        receiveNeighbor.coordinates(1) = (m_coordinateY + receiveNeighbor.direction(1) + m_nProcessorsY) % m_nProcessorsY;
        receiveNeighbor.coordinates(2) = (m_coordinateZ + receiveNeighbor.direction(2) + m_nProcessorsZ) % m_nProcessorsZ;

        sendNeighbor.rank = sendNeighbor.coordinates(0) * m_nProcessorsY * m_nProcessorsZ + sendNeighbor.coordinates(1) * m_nProcessorsZ + sendNeighbor.coordinates(2);

        receiveNeighbor.rank = sendNeighbor.coordinates(0) * m_nProcessorsY * m_nProcessorsZ + sendNeighbor.coordinates(1) * m_nProcessorsZ + sendNeighbor.coordinates(2);

        if(sendNeighbor.rank == world.rank()) {
            // If we're supposed to communicate with ourselves, we may
            // just push this information and skip the list of cells.
            m_sendNeighbors.push_back(sendNeighbor);
            m_receiveNeighbors.push_back(receiveNeighbor);
            continue;
        }

        // Find cells to transfer to this neighbor
        vector<MoleculeSystemCell*>& cellsToSend = sendNeighbor.cells;
        vector<MoleculeSystemCell*>& cellsToReceive = receiveNeighbor.cells;
        if(abs(direction(0)) == 1) { // x-directional faces
            int iSend;
            int iReceive;
            if(direction(0) == 1) {
                iSend = m_cellRangeX.lastElement();
                iReceive = iSend + 1;
            } else {
                iSend = m_cellRangeX.firstElement();
                iReceive = iSend - 1;
            }
            for(int jSend = m_cellRangeY.firstElement(); jSend <= m_cellRangeY.lastElement(); jSend++) {
                for(int kSend = m_cellRangeZ.firstElement(); kSend <= m_cellRangeZ.lastElement(); kSend++) {
                    int jReceive = jSend;
                    int kReceive = kSend;

                    int iSendLocal = (iSend + m_moleculeSystem->nCells()(0)) % m_moleculeSystem->nCells()(0);
                    int iReceiveLocal = (iReceive + m_moleculeSystem->nCells()(0)) % m_moleculeSystem->nCells()(0);
                    int jSendLocal = (jSend + m_moleculeSystem->nCells()(1)) % m_moleculeSystem->nCells()(1);
                    int jReceiveLocal = (jReceive + m_moleculeSystem->nCells()(1)) % m_moleculeSystem->nCells()(1);
                    int kSendLocal = (kSend + m_moleculeSystem->nCells()(2)) % m_moleculeSystem->nCells()(2);
                    int kReceiveLocal = (kReceive + m_moleculeSystem->nCells()(2)) % m_moleculeSystem->nCells()(2);
                    cellsToSend.push_back(m_moleculeSystem->cell(iSendLocal, jSendLocal, kSendLocal));
                    MoleculeSystemCell* cellToReceive = m_moleculeSystem->cell(iReceiveLocal, jReceiveLocal, kReceiveLocal);
                    cellsToReceive.push_back(cellToReceive);
                    if(std::find(m_localAndGhostCells.begin(), m_localAndGhostCells.end(), cellToReceive) == m_localAndGhostCells.end()) {
                        m_localAndGhostCells.push_back(cellToReceive);
                    }
                }
            }
        } else if(abs(direction(1)) == 1) { // y-directional faces
            int jSend;
            int jReceive;
            if(direction(1) == 1) {
                jSend = m_cellRangeY.lastElement();
                jReceive = jSend + 1;
            } else {
                jSend = m_cellRangeY.firstElement();
                jReceive = jSend - 1;
            }
            for(int iSend = m_cellRangeX.firstElement() - 1; iSend <= m_cellRangeX.lastElement() + 1; iSend++) {
                for(int kSend = m_cellRangeZ.firstElement(); kSend <= m_cellRangeZ.lastElement(); kSend++) {
                    int iReceive = iSend;
                    int kReceive = kSend;
                    int iSendLocal = (iSend + m_moleculeSystem->nCells()(0)) % m_moleculeSystem->nCells()(0);
                    int iReceiveLocal = (iReceive + m_moleculeSystem->nCells()(0)) % m_moleculeSystem->nCells()(0);
                    int jSendLocal = (jSend + m_moleculeSystem->nCells()(1)) % m_moleculeSystem->nCells()(1);
                    int jReceiveLocal = (jReceive + m_moleculeSystem->nCells()(1)) % m_moleculeSystem->nCells()(1);
                    int kSendLocal = (kSend + m_moleculeSystem->nCells()(2)) % m_moleculeSystem->nCells()(2);
                    int kReceiveLocal = (kReceive + m_moleculeSystem->nCells()(2)) % m_moleculeSystem->nCells()(2);
                    cellsToSend.push_back(m_moleculeSystem->cell(iSendLocal, jSendLocal, kSendLocal));
                    MoleculeSystemCell* cellToReceive = m_moleculeSystem->cell(iReceiveLocal, jReceiveLocal, kReceiveLocal);
                    cellsToReceive.push_back(cellToReceive);
                    if(std::find(m_localAndGhostCells.begin(), m_localAndGhostCells.end(), cellToReceive) == m_localAndGhostCells.end()) {
                        m_localAndGhostCells.push_back(cellToReceive);
                    }
                }
            }
        } else if(abs(direction(2)) == 1) { // z-directional faces
            int kSend;
            int kReceive;
            if(direction(2) == 1) {
                kSend = m_cellRangeZ.lastElement();
                kReceive = kSend + 1;
            } else {
                kSend = m_cellRangeZ.firstElement();
                kReceive = kSend - 1;
            }
            for(int iSend = m_cellRangeX.firstElement() - 1; iSend <= m_cellRangeX.lastElement() + 1; iSend++) {
                for(int jSend = m_cellRangeY.firstElement() - 1; jSend <= m_cellRangeY.lastElement() + 1; jSend++) {
                    int iReceive = iSend;
                    int jReceive = jSend;
                    int iSendLocal = (iSend + m_moleculeSystem->nCells()(0)) % m_moleculeSystem->nCells()(0);
                    int iReceiveLocal = (iReceive + m_moleculeSystem->nCells()(0)) % m_moleculeSystem->nCells()(0);
                    int jSendLocal = (jSend + m_moleculeSystem->nCells()(1)) % m_moleculeSystem->nCells()(1);
                    int jReceiveLocal = (jReceive + m_moleculeSystem->nCells()(1)) % m_moleculeSystem->nCells()(1);
                    int kSendLocal = (kSend + m_moleculeSystem->nCells()(2)) % m_moleculeSystem->nCells()(2);
                    int kReceiveLocal = (kReceive + m_moleculeSystem->nCells()(2)) % m_moleculeSystem->nCells()(2);
                    cellsToSend.push_back(m_moleculeSystem->cell(iSendLocal, jSendLocal, kSendLocal));
                    MoleculeSystemCell* cellToReceive = m_moleculeSystem->cell(iReceiveLocal, jReceiveLocal, kReceiveLocal);
                    cellsToReceive.push_back(cellToReceive);
                    if(std::find(m_localAndGhostCells.begin(), m_localAndGhostCells.end(), cellToReceive) == m_localAndGhostCells.end()) {
                        m_localAndGhostCells.push_back(cellToReceive);
                    }
                }
            }
        }
        m_sendNeighbors.push_back(sendNeighbor);
        m_receiveNeighbors.push_back(receiveNeighbor);
    }
    // We need to swap the receiving neighbors in each direction (x, y or z) because of the send/receive order
    swap(m_receiveNeighbors[0], m_receiveNeighbors[1]);
    swap(m_receiveNeighbors[2], m_receiveNeighbors[3]);
    swap(m_receiveNeighbors[4], m_receiveNeighbors[5]);
}

int Processor::rank()
{
    return world.rank();
}

void Processor::receiveAtomsFromNeighbor(const ProcessorNeighbor& neighbor) {
    //    int atomsRemoved = 0;
    //    cout << "About to receive atoms" << endl;
    for(MoleculeSystemCell* cellToReceive : neighbor.cells) {
//        cout << "Cell to receive atoms: " << cellToReceive->indices() << " has room for: " << cellToReceive->atoms().size() << endl;
        //        cout << "Receiving from cell" << endl;
        vector<Atom*> atomsToReceive; // IMPORTANT: This list must be freed manually! (Should be done below in this scope)
        pureCommunicationTimer.restart();
        world.recv(neighbor.rank, 0, atomsToReceive); // IMPORTANT: This list must be freed manually! (Should be done below in this scope)
//        cout << "Rank " << world.rank() << " receives from " << neighbor.rank << ", atoms for cell " << cellToReceive->indices() << ", in total:" << atomsToReceive.size() << endl;
        m_pureCommunicationTime += pureCommunicationTimer.elapsed();

//        cout << "Received " << atomsToReceive.size() << " atoms" << endl;

        uint nAtomsAvailable = cellToReceive->atoms().size();
        // Allocate space for new atoms if not enough atoms are available on this processor
        if(atomsToReceive.size() > nAtomsAvailable) {
            vector<Atom*> locallyAllocatedAtoms;
            for(uint i = nAtomsAvailable; i < atomsToReceive.size(); i++) {
                Atom* localAtom = new Atom();
                locallyAllocatedAtoms.push_back(localAtom);
            }
            cellToReceive->addAtoms(locallyAllocatedAtoms);
        } else if(atomsToReceive.size() < nAtomsAvailable) {
            for(uint i = atomsToReceive.size(); i < nAtomsAvailable; i++) {
                Atom* localAtom = cellToReceive->atoms().at(i);
                delete localAtom;
            }
            cellToReceive->deleteAtoms(nAtomsAvailable - atomsToReceive.size());
        }
        // For each received atom, make a clone to the local atom
        for(uint i = 0; i < atomsToReceive.size(); i++) {
            Atom* atom = atomsToReceive.at(i);
            Atom* localAtom = cellToReceive->atoms().at(i);
            // TODO Add a check to see if the particle types are actually set
            localAtom->communicationClone(*atom, m_moleculeSystem->particleTypes());
        }

        //  IMPORTANT: Need to explicitly free the memory used by MPI::Boost when allocating a received list of atoms
        for(Atom* atom : atomsToReceive) {
            delete atom;
        }
        atomsToReceive.clear();
    }
}

void Processor::sendAtomsToNeighbor(const ProcessorNeighbor& neighbor) {
    for(MoleculeSystemCell* cellToSend : neighbor.cells) {
//        cout << "Rank " << world.rank() << " sends to " << neighbor.rank << ", atoms from cell " << cellToSend->indices() << ", in total: " << cellToSend->atoms().size() << endl;
        vector<Atom*> atomsToSend;
        atomsToSend.insert(atomsToSend.end(), cellToSend->atoms().begin(), cellToSend->atoms().end());
        pureCommunicationTimer.restart();
        world.send(neighbor.rank, 0, atomsToSend);
        m_pureCommunicationTime += pureCommunicationTimer.elapsed();
        atomsToSend.clear();
    }
}

void Processor::clearForcesInNeighborCells() {
    for(const ProcessorNeighbor& neighbor : m_sendNeighbors) {
        for(MoleculeSystemCell* cell : neighbor.cells) {
            for(Atom* atom : cell->atoms()) {
                atom->clearForcePotentialPressure();
            }
        }
    }
    for(const ProcessorNeighbor& neighbor : m_receiveNeighbors) {
        for(MoleculeSystemCell* cell : neighbor.cells) {
            for(Atom* atom : cell->atoms()) {
                atom->clearForcePotentialPressure();
            }
        }
    }
}

bool Processor::shouldSendFirst(const irowvec& direction) {
    return (abs(direction(0)) && !(m_coordinateX % 2)) || (abs(direction(1)) && !(m_coordinateY % 2)) || (abs(direction(2)) && !(m_coordinateZ % 2));
}

int Processor::nAtoms()
{
    int numberOfAtoms = 0;
    for(MoleculeSystemCell* cell : m_localCells) {
        numberOfAtoms += cell->atoms().size();
    }
    return numberOfAtoms;
}

int Processor::nProcessors()
{
    return world.size();
}

void Processor::communicateAtoms() {
    m_communicationTimer.restart();
    // Send atoms to neighbors
    for(uint i = 0; i < m_directions.size(); i++) {
        irowvec direction = m_directions.at(i);
        const ProcessorNeighbor& sendNeighbor = m_sendNeighbors.at(i);
        const ProcessorNeighbor& receiveNeighbor = m_receiveNeighbors.at(i);

        if(shouldSendFirst(direction)) { // is even in the direction of transfer, send first, then receive
            sendAtomsToNeighbor(sendNeighbor);
            receiveAtomsFromNeighbor(receiveNeighbor);
        } else { // receive first, then send
            receiveAtomsFromNeighbor(receiveNeighbor);
            sendAtomsToNeighbor(sendNeighbor);
        }
    }
    m_totalCommunicationTime += m_communicationTimer.elapsed();
}

ProcessorNeighbor::ProcessorNeighbor() :
    direction(zeros<irowvec>(3)),
    coordinates(zeros<irowvec>(3))
{

}
