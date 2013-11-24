#include "twoparticleforce.h"

TwoParticleForce::TwoParticleForce() :
    m_isNewtonsThirdLawEnabled(false),
    m_isCalculatePressureEnabled(true),
    m_isCalculatePotentialEnabled(true),
    m_zeroVector(Vector3(0,0,0))
{
}

void TwoParticleForce::calculateAndApplyForce(Atom *atom1, Atom *atom2)
{
    calculateAndApplyForce(atom1, atom2, m_zeroVector);
}
