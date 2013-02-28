#include "andersenthermostat.h"

#include <src/moleculesystem.h>
#include <src/atom.h>
#include <src/molecule.h>
#include <src/random.h>
#include <src/integrator/integrator.h>

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
    for(Atom* atom : m_moleculeSystem->atoms()) {
        double randomNumber = random->ran2();
        if(randomNumber < dt / tau) {
            rowvec velocity = randn<rowvec>(m_nDimensions);
            velocity *= sqrt(targetTemperature / atom->mass());
            atom->setVelocity(velocity);
        }
    }
}

void AndersenThermostat::setTargetTemperature(double targetTemperature) {
    m_targetTemperature = targetTemperature;
}

void AndersenThermostat::setCollisionTime(double collisionTime) {
    m_collisionTime = collisionTime;
}

