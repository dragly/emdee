#include <src/integrator/integrator.h>
#include <src/atom.h>
#include <src/moleculesystem.h>

Integrator::Integrator(MoleculeSystem *moleculeSystem) :
    m_moleculeSystem(moleculeSystem),
    m_timeStep(0.01)
{
}

void Integrator::setTimeStep(double timeStep)
{
    m_timeStep = timeStep;
}
