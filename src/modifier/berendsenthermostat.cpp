#include "berendsenthermostat.h"

#include <src/moleculesystem.h>
#include <src/atom.h>
#include <src/integrator/integrator.h>

BerendsenThermostat::BerendsenThermostat(MoleculeSystem* moleculeSystem) :
    Modifier(moleculeSystem),
    m_targetTemperature(1.0),
    m_relaxationTime(1.0)
{
}

void BerendsenThermostat::apply()
{
    double dt = m_moleculeSystem->integrator()->timeStep();
    double tau = m_relaxationTime;
    double currentTemperature = m_moleculeSystem->temperature();
    double gamma = sqrt(1 + dt / tau * (m_targetTemperature / currentTemperature - 1));
    for(Atom* atom : m_moleculeSystem->atoms()) {
        atom->setVelocity(gamma * atom->velocity());
    }
}

void BerendsenThermostat::setTargetTemperature(double targetTemperature) {
    m_targetTemperature = targetTemperature;
}

void BerendsenThermostat::setRelaxationTime(double relaxationTime) {
    m_relaxationTime = relaxationTime;
}
