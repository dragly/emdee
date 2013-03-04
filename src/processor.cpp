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
//    vector<Atom> myAtoms;
    vector<Vector3> myPositions;
    vector<Vector3> myVelocities;
    vector<Vector3> myForces;
    for(MoleculeSystemCell* cell : m_cells) {
        for(Atom* atom : m_atoms) {
            myPositions.push_back(atom->position());
            myVelocities.push_back(atom->velocity());
            myForces.push_back(atom->force());
        }
    }
//    for(int i = 0; i < 10; i++) {
//        Atom* atom = new Atom(AtomType::argon());
////        Vector3* atom = new Vector3(world.rank(), i, 3); // = new Atom(AtomType::argon());
//        atom->setPosition(Vector3(world.rank(),i,3));
//        myAtoms.push_back(atom);
//    }
//    vector<Vector3*> allAtoms;
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

    cout << "I now have " << allAtoms.size() << " atoms" << endl;
    for(Atom* atom : allAtoms) {
        const AtomType& type = atom->type();
        AtomType type2 = atom->type();
        cout << type.abbreviation << endl;
        cout << type2.abbreviation << endl;
    }
//    allNumbers.reserve(1);

//    cout << "Trying to send " << myAtoms.size() << " atoms" << endl;

//    mpi::all_gather(world, &myAtoms[0], myAtoms.size(), allAtoms);

//    m_moleculeSystem->deleteAtoms();
//    m_moleculeSystem->addAtoms(allAtoms);
    m_moleculeSystem->refreshCellContents();
//    cout << "Banana!" << endl;
//    mpi::broadcast(world, allNumbers, 0);

//    for(Atom* atom : allAtoms) {
//    for(Vector3* atom : allAtoms) {
//        cout << "Found atom " << atom->position() << endl;
//        cout << "Found atom " << atom.position() << endl;
//    }
}
