#include "andersenthermostat.h"

#include <moleculesystem.h>
#include <atom.h>
#include <random.h>
#include <integrator/integrator.h>
#include <processor.h>
#include <moleculesystemcell.h>

/*!
 * \class AndersenThermostat
 * \brief The AndersenThermostat class modifies the temperature of the system by
 * use of a Andersen thermostat.
 */

AndersenThermostat::AndersenThermostat(MoleculeSystem* moleculeSystem) :
    Modifier(moleculeSystem),
    m_targetTemperature(1.0),
    m_collisionTime(1.0),
    m_nDimensions(3)
{
    random = new Random();
}

void AndersenThermostat::apply()
{
    double dt = m_moleculeSystem->integrator()->timeStep();
    double tau = m_collisionTime;
    double targetTemperature = m_targetTemperature;
    for(MoleculeSystemCell* cell : m_moleculeSystem->localCells()) {
        for(Atom* atom : cell->atoms()) {
            if(atom->isPositionFixed()) {
                continue;
            }
            double randomNumber = random->ran2();
            if(randomNumber < dt / tau) {
                Vector3 velocity(randn<rowvec>(m_nDimensions));
                velocity *= sqrt(targetTemperature / atom->mass());
                atom->setVelocity(velocity);
            }
        }
    }
}

void AndersenThermostat::setTargetTemperature(double targetTemperature) {
    m_targetTemperature = targetTemperature;
}

void AndersenThermostat::setCollisionTime(double collisionTime) {
    m_collisionTime = collisionTime;
}

