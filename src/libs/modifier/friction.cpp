#include "friction.h"

#include "moleculesystem.h"
#include "moleculesystemcell.h"
#include "math/vector3.h"
#include "atom.h"

/*!
 * \class Friction
 * \brief The Friction class applies a friction force to all particles of the system.
 */

Friction::Friction(MoleculeSystem *moleculeSystem) :
    Modifier(moleculeSystem)
{
}

void Friction::apply() {
    for(MoleculeSystemCell* cell : m_moleculeSystem->localCells()) {
        for(Atom* atom : cell->atoms()) {
            atom->addForce(-atom->velocity() * 2.0);
        }
    }
}
