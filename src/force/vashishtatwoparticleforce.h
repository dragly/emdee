#ifndef VASHISHTATWOPARTICLEFORCE_H
#define VASHISHTATWOPARTICLEFORCE_H

#include <src/force/twoparticleforce.h>

class VashishtaTwoParticleForce : public TwoParticleForce
{
public:
    VashishtaTwoParticleForce();

    void calculateAndApplyForce(Atom *atom1, Atom *atom2);
    void calculateAndApplyForce(Atom *atom1, Atom *atom2, const Vector3 &atomOffset);
};

#endif // VASHISHTATWOPARTICLEFORCE_H
