#ifndef VASHISHTATHREEPARTICLEFORCE_H
#define VASHISHTATHREEPARTICLEFORCE_H

#include <src/force/threeparticleforce.h>

class VashishtaThreeParticleForce : public ThreeParticleForce
{
public:
    VashishtaThreeParticleForce();

    void calculateAndApplyForce(Atom *atom1, Atom *atom2, Atom *atom3);
    void calculateAndApplyForce(Atom *atom1, Atom *atom2, Atom *atom3, const Vector3 &atom2Offset, const Vector3 &atom3Offset);
};

#endif // VASHISHTATHREEPARTICLEFORCE_H
