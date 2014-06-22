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

using namespace std;

/*!
 * \class Processor
 * \brief The Processor class holds information about the current processor and
 * neighboring processors. Its implementation depends on whether the program is
 * compiled with MPI or not.
 */

Processor::Processor(MoleculeSystem *moleculeSystem) :
    m_moleculeSystem(moleculeSystem),
    m_nProcessors(0),
    m_nProcessorsX(0),
    m_nProcessorsY(0),
    m_nProcessorsZ(0),
    m_totalCommunicationTime(0),
    m_pureCommunicationTime(0)
{
    LOG(INFO) << "Processor loaded. I am rank 0 of 1 processor (without MPI).";
}

void Processor::setupProcessors()
{
}

int Processor::rank()
{
    return 0;
}

void Processor::clearForcesInNeighborCells() {

}

bool Processor::shouldSendFirst(const irowvec& direction) {
    return false;
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
    return 1;
}

void Processor::communicateAtoms() {

}

void Processor::communicateForces() {

}

const vector<MoleculeSystemCell *> &Processor::localCells() const {
    return m_moleculeSystem->globalCells();
}

const vector<MoleculeSystemCell *> &Processor::localAndGhostCells() const
{
    return m_moleculeSystem->globalCells();
}
