#include "friction.h"

#include "moleculesystem.h"
#include "moleculesystemcell.h"
#include "math/vector3.h"
#include "atom.h"

Friction::Friction(MoleculeSystem *moleculeSystem) :
    Modifier(moleculeSystem)
{
}

void Friction::apply() {
    for(MoleculeSystemCell* cell : m_moleculeSystem->localCells()) {
        for(Atom* atom : cell->atoms()) {
            atom->addForce(-atom->velocity() * 0.1);
        }
    }
}
