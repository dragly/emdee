#include <integrator/integrator.h>

/*!
 * \class Integrator
 * \brief The Integrator class is the base class for all time integrator classes.
 */

Integrator::Integrator(MoleculeSystem *moleculeSystem) :
    m_moleculeSystem(moleculeSystem),
    m_timeStep(0.01)
{
}

void Integrator::setTimeStep(double timeStep)
{
    m_timeStep = timeStep;
}
