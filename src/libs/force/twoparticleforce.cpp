#include "twoparticleforce.h"

TwoParticleForce::TwoParticleForce() :
    m_isNewtonsThirdLawEnabled(false),
    m_isCalculatePressureEnabled(true),
    m_isCalculatePotentialEnabled(true),
    m_cutoffRadius(0.0)
{
}

void TwoParticleForce::calculateAndApplyForce(Atom *atom1, Atom *atom2)
{
    calculateAndApplyForce(atom1, atom2, Vector3::zeroVector());
}

double TwoParticleForce::cutoffRadius() const
{
    return m_cutoffRadius;
}

void TwoParticleForce::setCutoffRadius(double cutoffRadius)
{
    m_cutoffRadius = cutoffRadius;
}

