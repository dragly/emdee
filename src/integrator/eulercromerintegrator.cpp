#include <src/integrator/eulercromerintegrator.h>
#include <src/moleculesystem.h>
#include <src/molecule.h>


EulerCromerIntegrator::EulerCromerIntegrator(MoleculeSystem *moleculeSystem) :
    Integrator(moleculeSystem)
{
}

void EulerCromerIntegrator::initialize() {
    m_moleculeSystem->updateForces();
}

void EulerCromerIntegrator::stepForward() {
    double dt = m_timeStep;
    for(uint i = 0; i < m_moleculeSystem->atoms().size(); i++) {
        Atom *atom = m_moleculeSystem->atoms().at(i);
        rowvec velocity = atom->velocity();
        rowvec position = atom->position();
        velocity += atom->force() / (atom->mass()) * dt;
        atom->setVelocity(velocity);
        position += velocity * dt;
        atom->setPosition(position);
    }

    m_moleculeSystem->updateForces();
}
