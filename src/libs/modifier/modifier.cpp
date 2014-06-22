#include "modifier.h"

/*!
 * \class Modifier
 * \brief The Modifier class is the base class for all MoleculeSystem modifiers.
 */


Modifier::Modifier(MoleculeSystem *moleculeSystem) :
    m_moleculeSystem(moleculeSystem)
{
}
