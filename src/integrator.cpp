#include "integrator.h"
#include "atom.h"
#include "molecule.h"
#include "moleculesystem.h"

Integrator::Integrator(MoleculeSystem *moleculeSystem) :
    m_moleculeSystem(moleculeSystem),
    m_timeStep(0.01)
{
}

void Integrator::initialize() {
    m_moleculeSystem->updateForces();
}

void Integrator::setTimeStep(double timeStep)
{
    m_timeStep = timeStep;
}

void Integrator::stepForward() {
    double dt = m_timeStep;
    for(uint i = 0; i < m_moleculeSystem->molecules().size(); i++) {
        Molecule *molecule = m_moleculeSystem->molecules().at(i);
        rowvec velocity = molecule->velocity();
        rowvec position = molecule->position();
        velocity += molecule->force() / (2*molecule->mass()) * dt;
        molecule->setVelocity(velocity);
        position += velocity * dt;
        molecule->setPosition(position);
    }

    m_moleculeSystem->updateForces();

    for(uint i = 0; i < m_moleculeSystem->molecules().size(); i++) {
        Molecule *molecule = m_moleculeSystem->molecules().at(i);
        rowvec velocity = molecule->velocity();
        velocity += molecule->force() / (2*molecule->mass()) * dt;
        molecule->setVelocity(velocity);
    }
}
