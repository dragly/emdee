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
    m_nProcessorsZ(0)
{
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
    m_rangeX = Range(m_coordinateX, m_nProcessorsX, m_moleculeSystem->nCells()(0));
    m_rangeY = Range(m_coordinateY, m_nProcessorsY, m_moleculeSystem->nCells()(1));
    m_rangeZ = Range(m_coordinateZ, m_nProcessorsZ, m_moleculeSystem->nCells()(2));

    cout << "X Range: " << m_rangeX << endl;
    cout << "Y Range: " << m_rangeY << endl;
    cout << "Z Range: " << m_rangeZ << endl;

    for(int i = m_rangeX.firstElement(); i <= m_rangeX.lastElement(); i++) {
        for(int j = m_rangeY.firstElement(); j <= m_rangeY.lastElement(); j++) {
            for(int k = m_rangeZ.firstElement(); k <= m_rangeZ.lastElement(); k++) {
                MoleculeSystemCell* cell = m_moleculeSystem->cell(i,j,k);
                m_cells.push_back(cell);
            }
        }
    }

    // Find neighbors
}

int Processor::rank()
{
    return world.rank();
}

void Processor::communicateAtoms() {
    vector<Atom*> myAtoms;
    for(MoleculeSystemCell* cell : m_cells) {
        myAtoms.insert(myAtoms.end(), cell->atoms().begin(), cell->atoms().end());
    }
    vector<Atom*> allAtoms;
    if(world.rank() == 0) {
        allAtoms.insert(allAtoms.end(), myAtoms.begin(), myAtoms.end());

        for(int iRank = 1; iRank < world.size(); iRank++) {
            vector<Atom*> receivedAtoms;
            world.recv(iRank, 0, receivedAtoms);
            allAtoms.insert(allAtoms.begin(), receivedAtoms.begin(), receivedAtoms.end());
        }
        cout << "Received " << allAtoms.size() << " atoms" << endl;
        for(int iRank = 1; iRank < world.size(); iRank++) {
            world.send(iRank, 1, allAtoms);
        }
    } else {
        world.send(0, 0, myAtoms);
        world.recv(0, 1, allAtoms);
    }

    vector<Atom*> myLocalAtoms;

    cout << "I now have " << allAtoms.size() << " atoms" << endl;

    // Allocate new atoms because we don't trust Boost::MPI to do it for us
    for(Atom* atom : allAtoms) {
        Atom* localAtom = new Atom(AtomType::argon());
        localAtom->setPosition(atom->position());
        localAtom->setVelocity(atom->velocity());
        localAtom->setForce(atom->force());

        myLocalAtoms.push_back(localAtom);
    }

    m_moleculeSystem->deleteAtoms();
    m_moleculeSystem->addAtoms(myLocalAtoms);
    m_moleculeSystem->refreshCellContents();
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
