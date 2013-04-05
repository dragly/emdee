#ifndef SINGLEPARTICLEFORCE_H
#define SINGLEPARTICLEFORCE_H

class Atom;

class SingleParticleForce
{
public:
    SingleParticleForce();

    virtual void apply(Atom* atom) = 0;
};

#endif // SINGLEPARTICLEFORCE_H
