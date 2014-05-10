#ifndef KOHENTHREEPARTICLEFORCE_H
#define KOHENTHREEPARTICLEFORCE_H

#include "force/threeparticleforce.h"

/**
 * @brief The KohenThreeParticleForce class is not yet implemented and serves only as a benchmark test case.
 */
class KohenThreeParticleForce : public ThreeParticleForce
{
public:
    KohenThreeParticleForce();

    virtual void calculateAndApplyForce(Atom *atom1, Atom *atom2, Atom *atom3, const Vector3 &atom2Offset, const Vector3 &atom3Offset);
    virtual void calculateAndApplyForce(Atom *atom1, Atom *atom2, Atom *atom3);
private:
    Vector3 rij;
    Vector3 rik;
};

#endif // KOHENTHREEPARTICLEFORCE_H
