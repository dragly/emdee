#include <src/integrator/velocityverletintegrator.h>
#include <src/moleculesystem.h>
#include <src/atom.h>
#include <src/moleculesystemcell.h>
#include <src/processor.h>

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
    for(MoleculeSystemCell* cell : m_moleculeSystem->processor()->cells()) {
        for(Atom* atom : cell->atoms()) {
            if(atom->isPositionFixed()) {
                continue;
            }
            Vector3 velocity = atom->velocity();
            Vector3 position = atom->position();
            velocity += atom->force() / (2*atom->mass()) * dt;
            atom->setVelocity(velocity);
            position += velocity * dt;
            atom->setPosition(position);
        }
    }
    m_moleculeSystem->updateForces();

    for(MoleculeSystemCell* cell : m_moleculeSystem->processor()->cells()) {
        for(Atom* atom : cell->atoms()) {
            if(atom->isPositionFixed()) {
                continue;
            }
            Vector3 velocity = atom->velocity();
            velocity += atom->force() / (2*atom->mass()) * dt;
            atom->setVelocity(velocity);
        }
    }
}
