#include "berendsenthermostat.h"

#include <moleculesystem.h>
#include <atom.h>
#include <integrator/integrator.h>
#include <processor.h>
#include <moleculesystemcell.h>

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
    double argument = 1 + dt / tau * (m_targetTemperature / currentTemperature - 1);
    if(argument < 0) {
        return;
    }
    double gamma = sqrt(argument);
    for(MoleculeSystemCell* cell : m_moleculeSystem->localCells()) {
        for(Atom* atom : cell->atoms()) {
            if(atom->isPositionFixed()) {
                continue;
            }
            atom->setVelocity(gamma * atom->velocity());
        }
    }
}

void BerendsenThermostat::setTargetTemperature(double targetTemperature) {
    m_targetTemperature = targetTemperature;
}

void BerendsenThermostat::setRelaxationTime(double relaxationTime) {
    m_relaxationTime = relaxationTime;
}
