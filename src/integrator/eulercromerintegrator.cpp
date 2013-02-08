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
    for(uint i = 0; i < m_moleculeSystem->molecules().size(); i++) {
        Molecule *molecule = m_moleculeSystem->molecules().at(i);
        rowvec velocity = molecule->velocity();
        rowvec position = molecule->position();
        velocity += molecule->force() / (molecule->mass()) * dt;
        molecule->setVelocity(velocity);
        position += velocity * dt;
        molecule->setPosition(position);
    }

    m_moleculeSystem->updateForces();
}
