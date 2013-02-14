#include "berendsenthermostat.h"

#include <src/moleculesystem.h>
#include <src/atom.h>
#include <src/molecule.h>
#include <src/integrator/integrator.h>

BerendsenThermostat::BerendsenThermostat(MoleculeSystem* moleculeSystem) :
    Modifier(moleculeSystem)
{
}

void BerendsenThermostat::apply()
{
    double dt = m_moleculeSystem->integrator()->timeStep();
    double tau = 0.01;
    double targetTemperature = 0.5;
    double currentTemperature = m_moleculeSystem->temperature();
    double gamma = sqrt(1 + dt / tau * (targetTemperature / currentTemperature - 1));
    for(Molecule* molecule : m_moleculeSystem->molecules()) {
        molecule->setVelocity(gamma * molecule->velocity());
    }
}
