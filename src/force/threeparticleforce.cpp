#include "threeparticleforce.h"

ThreeParticleForce::ThreeParticleForce()
{
}

void ThreeParticleForce::calculateAndApplyForce(Atom *atom1, Atom *atom2, Atom *atom3)
{
    calculateAndApplyForce(atom1, atom2, atom3, m_zeroVector, m_zeroVector);
}
