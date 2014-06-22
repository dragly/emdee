#include "berendsenthermostat.h"

#include <moleculesystem.h>
#include <atom.h>
#include <integrator/integrator.h>
#include <processor.h>
#include <moleculesystemcell.h>

/*!
 * \class BerendsenThermostat
 * \brief The BerendsenThermostat class modifies the temperature of the system by
 * use of a Berendsen thermostat.
 */

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
    double gamma = dt / tau;
    double argument = 1 + gamma * (m_targetTemperature / currentTemperature - 1);
    if(argument < 0) {
        return;
    }
    double beta = sqrt(argument);
    for(MoleculeSystemCell* cell : m_moleculeSystem->localCells()) {
        for(Atom* atom : cell->atoms()) {
            if(atom->isPositionFixed()) {
                continue;
            }
            atom->setVelocity(beta * atom->velocity());
        }
    }
}

double BerendsenThermostat::targetTemperature() const
{
    return m_targetTemperature;
}

void BerendsenThermostat::setTargetTemperature(double targetTemperature)
{
    m_targetTemperature = targetTemperature;
}
double BerendsenThermostat::relaxationTime() const
{
    return m_relaxationTime;
}

void BerendsenThermostat::setRelaxationTime(double relaxationTime)
{
    m_relaxationTime = relaxationTime;
}



