#ifndef BERENDSENTHERMOSTAT_H
#define BERENDSENTHERMOSTAT_H

#include <modifier/modifier.h>

class BerendsenThermostat : public Modifier
{
public:
    BerendsenThermostat(MoleculeSystem *moleculeSystem);

    void apply();

    void setTargetTemperature(double targetTemperature);
    void setRelaxationTime(double relaxationTime);
protected:
    double m_targetTemperature;
    double m_relaxationTime;
};


#endif // BERENDSENTHERMOSTAT_H
