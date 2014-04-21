#ifndef FRICTION_H
#define FRICTION_H

#include "modifier/modifier.h"

class Friction : public Modifier
{
public:
    Friction(MoleculeSystem *moleculeSystem);

    void apply();
};

#endif // FRICTION_H
