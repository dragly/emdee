#ifndef ANDERSENTHERMOSTAT_H
#define ANDERSENTHERMOSTAT_H

#include <modifier/modifier.h>

class Random;

class AndersenThermostat : public Modifier
{
public:
    AndersenThermostat(MoleculeSystem *moleculeSystem);

    void apply();

    void setTargetTemperature(double targetTemperature);
    void setCollisionTime(double collisionTime);
protected:
    double m_targetTemperature;
    double m_collisionTime;
    Random *random;

    int m_nDimensions;
};

#endif // ANDERSENTHERMOSTAT_H
