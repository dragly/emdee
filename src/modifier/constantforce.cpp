#include "constantforce.h"

#include <src/moleculesystem.h>
#include <src/moleculesystemcell.h>
#include <src/atom.h>

ConstantForce::ConstantForce(MoleculeSystem* moleculeSystem) :
    Modifier(moleculeSystem)
{
}

void ConstantForce::apply()
{
    for(MoleculeSystemCell* cell : m_moleculeSystem->cells()) {
        for(Atom* atom : cell->atoms()) {
            if(atom->isPositionFixed()) {
                continue;
            }
            atom->addForce(m_forceVector);
        }
    }
}
