#ifndef CONSTANTFORCE_H
#define CONSTANTFORCE_H

#include <src/math/vector3.h>
#include <src/force/singleparticleforce.h>

class ConstantForce : public SingleParticleForce
{
public:
    ConstantForce();

    void apply(Atom *atom);
    void setForceVector(Vector3 forceVector);
private:
    Vector3 m_forceVector;
};

inline void ConstantForce::setForceVector(Vector3 forceVector) {
    m_forceVector = forceVector;
}

#endif // CONSTANTFORCE_H
