#include <src/integrator/eulercromerintegrator.h>
#include <src/moleculesystem.h>
#include <src/atom.h>
#include <src/processor.h>
#include <src/moleculesystemcell.h>


EulerCromerIntegrator::EulerCromerIntegrator(MoleculeSystem *moleculeSystem) :
    Integrator(moleculeSystem)
{
}

void EulerCromerIntegrator::initialize() {
    m_moleculeSystem->updateForces();
}

void EulerCromerIntegrator::stepForward() {
    double dt = m_timeStep;
    for(MoleculeSystemCell* cell : m_moleculeSystem->allCells()) {
        for(Atom* atom : cell->atoms()) {
            Vector3 velocity = atom->velocity();
            Vector3 position = atom->position();
            velocity += atom->force() / (atom->mass()) * dt;
            atom->setVelocity(velocity);
            position += velocity * dt;
            atom->setPosition(position);
        }
    }

    m_moleculeSystem->updateForces();
}
