#ifndef BERENDSENTHERMOSTAT_H
#define BERENDSENTHERMOSTAT_H

#include <src/modifier/modifier.h>

class BerendsenThermostat : public Modifier
{
public:
    BerendsenThermostat(MoleculeSystem *moleculeSystem);

    void apply();
};

#endif // BERENDSENTHERMOSTAT_H
