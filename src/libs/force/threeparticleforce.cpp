#include "threeparticleforce.h"

ThreeParticleForce::ThreeParticleForce() :
    m_isCalculatePotentialEnabled(false),
    m_isCalculatePressureEnabled(false),
    m_isNewtonsThirdLawEnabled(false)
{
}

void ThreeParticleForce::calculateAndApplyForce(Atom *atom1, Atom *atom2, Atom *atom3)
{
    calculateAndApplyForce(atom1, atom2, atom3, Vector3::zeroVector(), Vector3::zeroVector());
}

double ThreeParticleForce::cutoffRadius() const
{
    return m_cutoffRadius;
}

void ThreeParticleForce::setCutoffRadius(double cutoffRadius)
{
    m_cutoffRadius = cutoffRadius;
}
