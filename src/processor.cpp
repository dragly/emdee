#include "processor.h"
#include <src/moleculesystem.h>
#include <src/moleculesystemcell.h>
#include <src/range.h>
#include <src/atom.h>
#include <src/math/vector3.h>

#include <iostream>
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
    irowvec left = {-1, 0, 0};
    directions.push_back(left);
    irowvec right = {1, 0, 0};
    directions.push_back(right);
    irowvec front = {0, 1, 0};
    directions.push_back(front);
    irowvec back = {0, -1, 0};
    directions.push_back(back);
    irowvec down = {0, 0, -1};
    directions.push_back(down);
    irowvec up = {0, 0, 1};
    directions.push_back(up);

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
    m_nProcessors = world.size();
    cout << "Setting up " << m_nProcessors << " processors" << endl;
    m_nProcessorsZ = int( pow(m_nProcessors, 1./3.) + 1e-9);
    m_nProcessorsY = int( sqrt((double)m_nProcessors / (double)m_nProcessorsZ) + 1e-9);
    m_nProcessorsX = int( (double)m_nProcessors / ((double)m_nProcessorsZ * (double) m_nProcessorsY) + 1e-9);
    cout << "Dividing processor space into " << m_nProcessorsX << " * " << m_nProcessorsY << " * " << m_nProcessorsZ << " processors.";
    cout << "This gives " << m_moleculeSystem->nCells()(0) / (double)m_nProcessorsX << " * " << m_moleculeSystem->nCells()(1) / (double)m_nProcessorsY << " * " << m_moleculeSystem->nCells()(2) / (double)m_nProcessorsZ << " cells per processor" << endl;

    // Find own processor coordinates
    m_coordinateX = world.rank() / (m_nProcessorsZ * m_nProcessorsY);
    m_coordinateY = (world.rank() / m_nProcessorsZ) % m_nProcessorsY;
    m_coordinateZ = world.rank() % m_nProcessorsZ;

    cout << "My coordinate is " << m_coordinateX << " " << m_coordinateY << " " << m_coordinateZ << endl;

    // Find my cells
    m_cellRangeX = Range(m_coordinateX, m_nProcessorsX, m_moleculeSystem->nCells()(0));
    m_cellRangeY = Range(m_coordinateY, m_nProcessorsY, m_moleculeSystem->nCells()(1));
    m_cellRangeZ = Range(m_coordinateZ, m_nProcessorsZ, m_moleculeSystem->nCells()(2));

    cout << "X Range: " << m_cellRangeX << endl;
    cout << "Y Range: " << m_cellRangeY << endl;
    cout << "Z Range: " << m_cellRangeZ << endl;

    for(int i = m_cellRangeX.firstElement(); i <= m_cellRangeX.lastElement(); i++) {
        for(int j = m_cellRangeY.firstElement(); j <= m_cellRangeY.lastElement(); j++) {
            for(int k = m_cellRangeZ.firstElement(); k <= m_cellRangeZ.lastElement(); k++) {
                //                cout << "I have cells " << i << " " << j << " " << k << endl;
                MoleculeSystemCell* cell = m_moleculeSystem->cell(i,j,k);
                if(m_cellRangeX.firstElement() == i || m_cellRangeY.firstElement() == j || m_cellRangeZ.firstElement() == k ||
                        m_cellRangeX.lastElement() == i || m_cellRangeY.lastElement() == j || m_cellRangeZ.lastElement() == k) {
                    cell->setOnProcessorEdge(true);
                }
                m_cells.push_back(cell);
            }
        }
    }

    // Find neighbors
    for(const irowvec& direction : directions) {

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
            cout << "No need to communicate with self..." << endl;
            sendNeighbors.push_back(sendNeighbor);
            receiveNeighbors.push_back(receiveNeighbor);
            continue;
        }

//        cout << "Neighbor in direction " << sendNeighbor.direction;

//        cout << "Neighbor at " << sendNeighbor.coordinates;

//        cout << "Has ID " << sendNeighbor.rank << endl;

        // Find cells to transfer to this neighbor
        vector<MoleculeSystemCell*>& cellsToSend = sendNeighbor.cells;
        vector<MoleculeSystemCell*>& cellsToReceive = receiveNeighbor.cells;
        if(abs(direction(0)) == 1) { // x-directional faces
//            cout << "X-directional send/recv";
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
                    //                            cout << "Wanna transfer cell " << iLocal << "  " << jLocal << " " << kLocal << endl;
                    cellsToSend.push_back(m_moleculeSystem->cell(iSendLocal, jSendLocal, kSendLocal));
                    cellsToReceive.push_back(m_moleculeSystem->cell(iReceiveLocal, jReceiveLocal, kReceiveLocal));
//                    cout << "send" << m_moleculeSystem->cell(iSendLocal, jSendLocal, kSendLocal)->indices();
//                    cout << "recv" << m_moleculeSystem->cell(iReceiveLocal, jReceiveLocal, kReceiveLocal)->indices();
                }
            }
        } else if(abs(direction(1)) == 1) { // y-directional faces
//            cout << "Y-directional send/recv";
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
                    //                            cout << "Wanna transfer cell " << iLocal << "  " << jLocal << " " << kLocal << endl;
                    cellsToSend.push_back(m_moleculeSystem->cell(iSendLocal, jSendLocal, kSendLocal));
                    cellsToReceive.push_back(m_moleculeSystem->cell(iReceiveLocal, jReceiveLocal, kReceiveLocal));
//                    cout << "send" << m_moleculeSystem->cell(iSendLocal, jSendLocal, kSendLocal)->indices();
//                    cout << "recv" << m_moleculeSystem->cell(iReceiveLocal, jReceiveLocal, kReceiveLocal)->indices();
                }
            }
        } else if(abs(direction(2)) == 1) { // z-directional faces
//            cout << "Z-directional send/recv";
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
                    //                            cout << "Wanna transfer cell " << iLocal << "  " << jLocal << " " << kLocal << endl;
                    cellsToSend.push_back(m_moleculeSystem->cell(iSendLocal, jSendLocal, kSendLocal));
                    cellsToReceive.push_back(m_moleculeSystem->cell(iReceiveLocal, jReceiveLocal, kReceiveLocal));
//                    cout << "send" << m_moleculeSystem->cell(iSendLocal, jSendLocal, kSendLocal)->indices();
//                    cout << "recv" << m_moleculeSystem->cell(iReceiveLocal, jReceiveLocal, kReceiveLocal)->indices();
                }
            }
        }
//        cout << "Wanna send " << cellsToSend.size() << " cells" << endl;
//        cout << "Wanna receive " << cellsToReceive.size() << " cells" << endl;

        cout << endl;
        sendNeighbors.push_back(sendNeighbor);
        receiveNeighbors.push_back(receiveNeighbor);
    }


//    cout << "Setting up force neighbors" << endl;
//    vector<Range> cellRanges;
//    cellRanges.push_back(m_cellRangeX);
//    cellRanges.push_back(m_cellRangeY);
//    cellRanges.push_back(m_cellRangeZ);
//    // Find neighbors
//    for(const irowvec& direction : forceDirections) {
//        for(int isSend = 0; isSend < 2; isSend++) {
//            ProcessorNeighbor neighbor;
//            if(isSend == 1) {
//                neighbor.direction = direction;
//            } else {
//                neighbor.direction = -direction;
//            }

//            int i = neighbor.direction(0);
//            int j = neighbor.direction(1);
//            int k = neighbor.direction(2);

//            neighbor.coordinates(0) = (m_coordinateX + i + m_nProcessorsX) % m_nProcessorsX;
//            neighbor.coordinates(1) = (m_coordinateY + j + m_nProcessorsY) % m_nProcessorsY;
//            neighbor.coordinates(2) = (m_coordinateZ + k + m_nProcessorsZ) % m_nProcessorsZ;

//            neighbor.rank = neighbor.coordinates(0) * m_nProcessorsY * m_nProcessorsZ + neighbor.coordinates(1) * m_nProcessorsZ + neighbor.coordinates(2);

//            if(neighbor.rank == world.rank()) {
//                cout << "No need to communicate with self..." << endl;
//                if(isSend == 1) {
//                    forceSendNeighbors.push_back(neighbor);
//                } else {
//                    forceReceiveNeighbors.push_back(neighbor);
//                }
//                continue;
//            }

//            cout << "In direction ";
//            neighbor.direction.print();
//            cout << " there is a neighbor at " << neighbor.coordinates;

//            cout << "Neighbor has ID " << neighbor.rank << endl;

//            // Find cells to transfer to this neighbor
//            vector<MoleculeSystemCell*>& cells = neighbor.cells;

//            if(isSend == 1) {
//                cout << "Sending" << endl;
//            } else {
//                cout << "Receiving" << endl;
//            }
//            imat minMax = zeros<imat>(3,2);

//            for(int iDim = 0; iDim < 3; iDim++) {
//                if(abs(direction(iDim)) == 0) {
//                    minMax(iDim, 0) = cellRanges[iDim].firstElement();
//                    minMax(iDim, 1) = cellRanges[iDim].lastElement();
//                } else {
//                    if(isSend == 1) {
//                        if(direction(iDim) == 1) {
//                            minMax(iDim, 0) = cellRanges[iDim].lastElement();
//                            minMax(iDim, 1) = cellRanges[iDim].lastElement();
//                        } else {
//                            minMax(iDim, 0) = cellRanges[iDim].firstElement();
//                            minMax(iDim, 1) = cellRanges[iDim].firstElement();
//                        }
//                    } else {
//                        if(direction(iDim) == 1) {
//                            minMax(iDim, 0) = cellRanges[iDim].firstElement() - 1;
//                            minMax(iDim, 1) = cellRanges[iDim].firstElement() - 1;
//                        } else {
//                            minMax(iDim, 0) = cellRanges[iDim].lastElement() + 1;
//                            minMax(iDim, 1) = cellRanges[iDim].lastElement() + 1;
//                        }
//                    }
//                }
//            }

////            cout << "MIN MAX" << endl << minMax;
//            for(int i = minMax(0, 0); i <= minMax(0, 1); i++) {
//                for(int j = minMax(1, 0); j <= minMax(1, 1); j++) {
//                    for(int k = minMax(2, 0); k <= minMax(2, 1); k++) {
////                        cout << "Pushing cell " << m_moleculeSystem->cell(i,j,k)->indices();
//                        cells.push_back(m_moleculeSystem->cell(i,j,k));
//                    }
//                }
//            }

//            if(isSend == 1) {
////                cout << "Wanna send " << cells.size() << " cells" << endl;
//                forceSendNeighbors.push_back(neighbor);
//            } else {
////                cout << "Wanna receive " << cells.size() << " cells" << endl;
//                forceReceiveNeighbors.push_back(neighbor);
//            }
//            cout << "Done..." << endl;
//        }
//    }

}

int Processor::rank()
{
    return world.rank();
}

void Processor::receiveAtomsFromNeighbor(const ProcessorNeighbor& neighbor) {
    int atomsRemoved = 0;
    for(MoleculeSystemCell* cellToReceive : neighbor.cells) {
        vector<Atom*> atomsToReceive;
//        cout << "Removed " << cellToReceive->atoms().size() << " atoms from " << cellToReceive->indices();
        atomsRemoved += cellToReceive->atoms().size();
        cellToReceive->deleteAtomsFromCellAndSystem();
        pureCommunicationTimer.restart();
        world.recv(neighbor.rank, 0, atomsToReceive);
        m_pureCommunicationTime += pureCommunicationTimer.elapsed();
        vector<Atom*> locallyAllocatedAtoms;
        for(Atom* atom : atomsToReceive) {
            Atom* localAtom = new Atom(AtomType::argon());
//            localAtom->clone(*atom);
            localAtom->communicationClone(*atom);
            locallyAllocatedAtoms.push_back(localAtom);
        }
        m_moleculeSystem->addAtoms(locallyAllocatedAtoms);
        cellToReceive->addAtoms(locallyAllocatedAtoms);
    }
//    cout << "Removed a total of " << atomsRemoved << " atoms" << endl;
}

void Processor::sendAtomsToNeighbor(const ProcessorNeighbor& neighbor) {
    for(MoleculeSystemCell* cellToSend : neighbor.cells) {
        vector<Atom*> atomsToSend;
        atomsToSend.insert(atomsToSend.end(), cellToSend->atoms().begin(), cellToSend->atoms().end());
//        cout << "Sending " << cellToSend->atoms().size() << " atoms from " << cellToSend->indices();
        pureCommunicationTimer.restart();
        world.send(neighbor.rank, 0, atomsToSend);
        m_pureCommunicationTime += pureCommunicationTimer.elapsed();
    }
}

bool Processor::shouldSendFirst(const irowvec& direction) {
    return (abs(direction(0)) && !(m_coordinateX % 2)) || (abs(direction(1)) && !(m_coordinateY % 2)) || (abs(direction(2)) && !(m_coordinateZ % 2));
}

void Processor::communicateAtoms() {
    communicationTimer.restart();
    //    m_moleculeSystem->refreshCellContents();

    // Send atoms to neighbors
    for(uint i = 0; i < directions.size(); i++) {
        irowvec direction = directions.at(i);
        const ProcessorNeighbor& sendNeighbor = sendNeighbors.at(i);
        const ProcessorNeighbor& receiveNeighbor = receiveNeighbors.at(i);
//            cout << "Supposed to communicate with " << neighbor.rank << " in direction " << neighbor.direction << endl;

//            cout << "Started out with " << m_moleculeSystem->atoms().size() << " atoms in the system" << endl;

        if(shouldSendFirst(direction)) { // is even in the direction of transfer, send first, then receive
            sendAtomsToNeighbor(sendNeighbor);
            receiveAtomsFromNeighbor(receiveNeighbor);
        } else { // receive first, then send
            receiveAtomsFromNeighbor(receiveNeighbor);
            sendAtomsToNeighbor(sendNeighbor);
        }
//            cout << "Removed atoms from cell, now I have " << m_moleculeSystem->atoms().size() << " atoms in the system" << endl;

//            cout << "Currently have " << m_moleculeSystem->atoms().size() << " atoms in the system" << endl;

//            m_moleculeSystem->addAtomsToCorrectCells(locallyAllocatedAtoms);
    }
    cout << "Now I have " << m_moleculeSystem->atoms().size() << " atoms in the system" << endl;
    m_totalCommunicationTime += communicationTimer.elapsed();
    cout << "Total communication time so far: " << m_totalCommunicationTime << endl;
    cout << "Pure communication time so far: " << m_pureCommunicationTime << endl;
}


//void Processor::receiveForcesFromNeighbor(const ProcessorNeighbor& neighbor) {
//    //    int atomsRemoved = 0;
//    //    cout << "Expecting to receive " << neighbor.forceCells.size() << " cells" << endl;
//    for(MoleculeSystemCell* cellToReceive : neighbor.cells) {
//        vector<Atom*> atomsToReceive;
//        //        cout << "Removed " << cellToReceive->atoms().size() << " atoms from " << cellToReceive->indices();
//        //        atomsRemoved += cellToReceive->atoms().size();
//        //        cellToReceive->deleteAtomsFromCellAndSystem();
//        pureCommunicationTimer.restart();
//        world.recv(neighbor.rank, 0, atomsToReceive);
//        m_pureCommunicationTime += pureCommunicationTimer.elapsed();
//        //        vector<Atom*> locallyAllocatedAtoms;
//        cout << "Receiving cell " << cellToReceive->indices() << " which should have " << cellToReceive->atoms().size() << " atoms" << "and I will be receiving " << atomsToReceive.size() << " atoms" << endl;
//        for(uint iAtom = 0; iAtom < cellToReceive->atoms().size(); iAtom++) {
//            //            Atom* localAtom = new Atom(AtomType::argon());
//            //            localAtom->clone(*atom);
//            //            localAtom->communicationClone(*atom);
//            Atom* localAtom = cellToReceive->atoms().at(iAtom);
//            Atom* receivedAtom = atomsToReceive.at(iAtom);
//            localAtom->addForce(receivedAtom->force());
//            //            locallyAllocatedAtoms.push_back(localAtom);
//        }
//        //        m_moleculeSystem->addAtoms(locallyAllocatedAtoms);
//        //        cellToReceive->addAtoms(locallyAllocatedAtoms);
//    }
//    //    cout << "Removed a total of " << atomsRemoved << " atoms" << endl;
//    //    cout << "Cells received" << endl;
//}

//void Processor::sendForcesToNeighbor(const ProcessorNeighbor& neighbor) {
//    //    cout << "Expecting to send " << neighbor.cellsToSendForces.size() << " cells" << endl;
//    for(MoleculeSystemCell* cellToSend : neighbor.cells) {
//        vector<Atom*> atomsToSend;
//        //        cout << "Sending cell cell " << cellToSend->indices() << endl;
//        atomsToSend.insert(atomsToSend.end(), cellToSend->atoms().begin(), cellToSend->atoms().end());
//        cout << "Sending  from cell " <<  cellToSend->indices() << " with " << cellToSend->atoms().size() << " atoms" << endl;
//        pureCommunicationTimer.restart();
//        world.send(neighbor.rank, 0, atomsToSend);
//        m_pureCommunicationTime += pureCommunicationTimer.elapsed();
//    }
//    //    cout << "Cells sent" << endl;
//}

//void Processor::communicateForces() {
//    cout << "Started force communication" << endl;
//    for(uint i = 0; i < forceDirections.size(); i++) {
//        cout << "Communication direction is " << forceDirections.at(i) << endl;
//        const ProcessorNeighbor& sendNeighbor = forceSendNeighbors.at(i);
//        const ProcessorNeighbor& receiveNeighbor = forceReceiveNeighbors.at(i);
//        cout << "Communicating with rank " << sendNeighbor.rank << " (I am " << world.rank() << ") " << endl;
//        if(sendNeighbor.rank == world.rank()) {
//            //            cout << "No need to communicate with self..." << endl;
//            continue;
//        }
//        if(shouldSendFirst(forceDirections.at(i))) {
//            cout << "Should send with " << sendNeighbor.rank << endl;
//            cout << "Should receive with " << receiveNeighbor.rank << endl;
//            sendForcesToNeighbor(sendNeighbor);
//            receiveForcesFromNeighbor(receiveNeighbor);
//        } else {
//            cout << "Should receive with " << receiveNeighbor.rank << endl;
//            cout << "Should send with " << sendNeighbor.rank << endl;
//            receiveForcesFromNeighbor(sendNeighbor);
//            sendForcesToNeighbor(receiveNeighbor);
//        }
//    }
//}


ProcessorNeighbor::ProcessorNeighbor() :
    direction(zeros<irowvec>(3)),
    coordinates(zeros<irowvec>(3))
{

}