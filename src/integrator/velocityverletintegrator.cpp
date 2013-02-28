#include <src/integrator/velocityverletintegrator.h>
#include <src/moleculesystem.h>
#include <src/atom.h>

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
    for(uint i = 0; i < m_moleculeSystem->atoms().size(); i++) {
        Atom *atom = m_moleculeSystem->atoms().at(i);
        rowvec velocity = atom->velocity();
        rowvec position = atom->position();
        velocity += atom->force() / (2*atom->mass()) * dt;
        atom->setVelocity(velocity);
        position += velocity * dt;
        atom->setPosition(position);
    }
    m_moleculeSystem->updateForces();

    for(uint i = 0; i < m_moleculeSystem->atoms().size(); i++) {
        Atom *atom = m_moleculeSystem->atoms().at(i);
        rowvec velocity = atom->velocity();
        velocity += atom->force() / (2*atom->mass()) * dt;
        atom->setVelocity(velocity);
    }
}
