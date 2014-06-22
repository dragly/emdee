#include "constantforce.h"

#include <atom.h>

/*!
 * \class ConstantForce
 * \brief The ConstantForce class defines a simple constant force that
 * is applied to all particles.
 */

ConstantForce::ConstantForce()
{
}

void ConstantForce::apply(Atom *atom)
{
    if(atom->isPositionFixed()) {
        return;
    }
    atom->addForce(m_forceVector);
}
