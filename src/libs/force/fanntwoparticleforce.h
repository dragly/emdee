#ifndef FANNTWOPARTICLEFORCE_H
#define FANNTWOPARTICLEFORCE_H

#include <force/twoparticleforce.h>

class FannTwoParticleForce : public TwoParticleForce
{
public:
    FannTwoParticleForce();

    // TwoParticleForce interface
public:
    virtual void calculateAndApplyForce(Atom *atom1, Atom *atom2);
    virtual void calculateAndApplyForce(Atom *atom1, Atom *atom2, const Vector3 &atom2Offset);
};

#endif // FANNTWOPARTICLEFORCE_H
