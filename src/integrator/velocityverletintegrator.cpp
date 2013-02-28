#include <src/integrator/velocityverletintegrator.h>
#include <src/moleculesystem.h>
#include <src/molecule.h>

VelocityVerletIntegrator::VelocityVerletIntegrator(MoleculeSystem *moleculeSystem) :
    Integrator(moleculeSystem)
{
}

void VelocityVerletIntegrator::initialize() {
    m_moleculeSystem->updateForces();
}

// TODO Build the timestep into the velocity to reduce the number of computations
void VelocityVerletIntegrator::stepForward() {
    double dt = m_timeStep;
    for(uint i = 0; i < m_moleculeSystem->molecules().size(); i++) {
        Atom *molecule = m_moleculeSystem->molecules().at(i);
        rowvec velocity = molecule->velocity();
        rowvec position = molecule->position();
        velocity += molecule->force() / (2*molecule->mass()) * dt;
        molecule->setVelocity(velocity);
        position += velocity * dt;
        molecule->setPosition(position);
    }
    m_moleculeSystem->updateForces();

    for(uint i = 0; i < m_moleculeSystem->molecules().size(); i++) {
        Atom *molecule = m_moleculeSystem->molecules().at(i);
        rowvec velocity = molecule->velocity();
        velocity += molecule->force() / (2*molecule->mass()) * dt;
        molecule->setVelocity(velocity);
    }
}
