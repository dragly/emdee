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
    m_totalCommunicationTime(0)
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
                m_cells.push_back(cell);
            }
        }
    }

    // Find neighbors
    for(const irowvec& direction : directions) {

        ProcessorNeighbor neighbor;

        neighbor.direction = direction;

        int i = direction(0);
        int j = direction(1);
        int k = direction(2);

        neighbor.coordinates(0) = (m_coordinateX + i + m_nProcessorsX) % m_nProcessorsX;
        neighbor.coordinates(1) = (m_coordinateY + j + m_nProcessorsY) % m_nProcessorsY;
        neighbor.coordinates(2) = (m_coordinateZ + k + m_nProcessorsZ) % m_nProcessorsZ;

        neighbor.rank = neighbor.coordinates(0) * m_nProcessorsY * m_nProcessorsZ + neighbor.coordinates(1) * m_nProcessorsZ + neighbor.coordinates(2);

        if(neighbor.rank == world.rank()) {
            cout << "No need to communicate with self..." << endl;
            continue;
        }

        cout << "Neighbor in direction " << neighbor.direction;

        cout << "Neighbor at " << neighbor.coordinates;

        cout << "Has ID " << neighbor.rank << endl;

        // Find cells to transfer to this neighbor
        vector<MoleculeSystemCell*>& cellsToSend = neighbor.cellsToSend;
        vector<MoleculeSystemCell*>& cellsToReceive = neighbor.cellsToReceive;
        if(abs(i) == 1) { // x-directional faces
//            cout << "X-directional send/recv";
            int iSend;
            int iReceive;
            if(i == 1) {
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
        }
        if(abs(j) == 1) { // y-directional faces
//            cout << "Y-directional send/recv";
            int jSend;
            int jReceive;
            if(j == 1) {
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
        }
        if(abs(k) == 1) { // z-directional faces
//            cout << "Z-directional send/recv";
            int kSend;
            int kReceive;
            if(k == 1) {
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
        cout << "Wanna send " << cellsToSend.size() << " cells" << endl;
        cout << "Wanna receive " << cellsToReceive.size() << " cells" << endl;

        cout << endl;
        processorNeighbors.push_back(neighbor);
    }
}

int Processor::rank()
{
    return world.rank();
}

void Processor::receiveAtomsFromNeighbor(const ProcessorNeighbor& neighbor) {
    int atomsRemoved = 0;
    for(MoleculeSystemCell* cellToReceive : neighbor.cellsToReceive) {
        vector<Atom*> atomsToReceive;
//        cout << "Removed " << cellToReceive->atoms().size() << " atoms from " << cellToReceive->indices();
        atomsRemoved += cellToReceive->atoms().size();
        cellToReceive->deleteAtomsFromCellAndSystem();
        world.recv(neighbor.rank, 0, atomsToReceive);
        vector<Atom*> locallyAllocatedAtoms;
        for(Atom* atom : atomsToReceive) {
            Atom* localAtom = new Atom(AtomType::argon());
            localAtom->clone(*atom);
            locallyAllocatedAtoms.push_back(localAtom);
        }
        m_moleculeSystem->addAtoms(locallyAllocatedAtoms);
        cellToReceive->addAtoms(locallyAllocatedAtoms);
    }
//    cout << "Removed a total of " << atomsRemoved << " atoms" << endl;
}

void Processor::sendAtomsToNeighbor(const ProcessorNeighbor& neighbor) {
    for(MoleculeSystemCell* cellToSend : neighbor.cellsToSend) {
        vector<Atom*> atomsToSend;
        atomsToSend.insert(atomsToSend.end(), cellToSend->atoms().begin(), cellToSend->atoms().end());
//        cout << "Sending " << cellToSend->atoms().size() << " atoms from " << cellToSend->indices();
        world.send(neighbor.rank, 0, atomsToSend);
    }
}

void Processor::communicateAtoms() {
    communicationTimer.restart();
    //    m_moleculeSystem->refreshCellContents();

    // Send atoms to neighbors
    for(const irowvec& direction : directions) {
        for(const ProcessorNeighbor& neighbor : processorNeighbors) {
            urowvec equal = (neighbor.direction == direction);
            if(min(equal) == 0) {
                continue;
            }
//            cout << "Supposed to communicate with " << neighbor.rank << " in direction " << neighbor.direction << endl;

//            cout << "Started out with " << m_moleculeSystem->atoms().size() << " atoms in the system" << endl;

            if((abs(neighbor.direction(0)) && !(m_coordinateX % 2)) || (abs(neighbor.direction(1)) && !(m_coordinateY % 2)) || (abs(neighbor.direction(2)) && !(m_coordinateZ % 2))) { // is even in the direction of transfer, send first, then receive
                sendAtomsToNeighbor(neighbor);
                receiveAtomsFromNeighbor(neighbor);
            } else { // receive first, then send
                receiveAtomsFromNeighbor(neighbor);
                sendAtomsToNeighbor(neighbor);
            }
//            cout << "Removed atoms from cell, now I have " << m_moleculeSystem->atoms().size() << " atoms in the system" << endl;

//            cout << "Currently have " << m_moleculeSystem->atoms().size() << " atoms in the system" << endl;

//            m_moleculeSystem->addAtomsToCorrectCells(locallyAllocatedAtoms);

        }
    }
    cout << "Now I have " << m_moleculeSystem->atoms().size() << " atoms in the system" << endl;
    m_totalCommunicationTime += communicationTimer.elapsed();
    cout << "Total communication time so far: " << m_totalCommunicationTime << endl;
//    exit(999);
    // Receive atoms from neighbors

    //    vector<Atom*> myAtoms;
    //    for(MoleculeSystemCell* cell : m_cells) {
    //        myAtoms.insert(myAtoms.end(), cell->atoms().begin(), cell->atoms().end());
    //    }
    //    vector<Atom*> allAtoms;
    //    if(world.rank() == 0) {
    //        allAtoms.insert(allAtoms.end(), myAtoms.begin(), myAtoms.end());

    //        for(int iRank = 1; iRank < world.size(); iRank++) {
    //            vector<Atom*> receivedAtoms;
    //            world.recv(iRank, 0, receivedAtoms);
    //            allAtoms.insert(allAtoms.begin(), receivedAtoms.begin(), receivedAtoms.end());
    //        }
    //        for(int iRank = 1; iRank < world.size(); iRank++) {
    //            world.send(iRank, 1, allAtoms);
    //        }
    //    } else {
    //        world.send(0, 0, myAtoms);
    //        world.recv(0, 1, allAtoms);
    //    }

    //    vector<Atom*> myLocalAtoms;

    //    // Allocate new atoms because we don't trust Boost::MPI to do it for us
    //    for(uint i = 0; i < allAtoms.size(); i++) {
    //        Atom* atom = allAtoms.at(i);
    //        Atom* localAtom = new Atom(AtomType::argon());
    //        localAtom->clone(*atom);
    //        myLocalAtoms.push_back(localAtom);
    //    }
    //    m_moleculeSystem->deleteAtoms();
    //    m_moleculeSystem->addAtoms(myLocalAtoms);


}


//vector<Vector3> myPositions;
//vector<Vector3> myVelocities;
//vector<Vector3> myForces;
//for(MoleculeSystemCell* cell : m_cells) {
//    for(Atom* atom : cell->atoms()) {
//        myPositions.push_back(atom->position());
//        myVelocities.push_back(atom->velocity());
//        myForces.push_back(atom->force());
//    }
//}

//vector<Vector3> allPositions;
//vector<Vector3> allVelocities;
//vector<Vector3> allForces;
//if(world.rank() == 0) {
//    allPositions.insert(allPositions.end(), myPositions.begin(), myPositions.end());
//    allVelocities.insert(allVelocities.end(), myVelocities.begin(), myVelocities.end());
//    allForces.insert(allForces.end(), myForces.begin(), myForces.end());

//    for(int iRank = 1; iRank < world.size(); iRank++) {
//        vector<Vector3> receivedPositions;
//        vector<Vector3> receivedVelocities;
//        vector<Vector3> receivedForces;
//        world.recv(iRank, 0, receivedPositions);
//        world.recv(iRank, 1, receivedVelocities);
//        world.recv(iRank, 2, receivedForces);
//        allPositions.insert(allPositions.begin(), receivedPositions.begin(), receivedPositions.end());
//        allVelocities.insert(allVelocities.begin(), receivedVelocities.begin(), receivedVelocities.end());
//        allForces.insert(allForces.begin(), receivedForces.begin(), receivedForces.end());
//    }
//    cout << "Received " << allPositions.size() << " positions" << endl;
//    for(int iRank = 1; iRank < world.size(); iRank++) {
//        world.send(iRank, 10, allPositions);
//        world.send(iRank, 11, allVelocities);
//        world.send(iRank, 12, allForces);
//    }
//} else {
//    world.send(0, 0, myPositions);
//    world.send(0, 1, myVelocities);
//    world.send(0, 2, myForces);
//    world.recv(0, 10, allPositions);
//    world.recv(0, 11, allVelocities);
//    world.recv(0, 12, allForces);
//}

//vector<Atom*> allAtoms;

//for(int i = 0; i < allPositions.size(); i++) {
//    Atom* atom = new Atom(AtomType::argon());
//    atom->setPosition(allPositions.at(i));
//    atom->setVelocity(allVelocities.at(i));
//    atom->setForce(allForces.at(i));
//}


ProcessorNeighbor::ProcessorNeighbor() :
    direction(zeros<irowvec>(3)),
    coordinates(zeros<irowvec>(3))
{

}
