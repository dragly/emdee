#ifndef CONSTANTFORCE_H
#define CONSTANTFORCE_H

#include <src/modifier/modifier.h>
#include <src/math/vector3.h>

class ConstantForce : public Modifier
{
public:
    ConstantForce(MoleculeSystem *moleculeSystem);
    void apply();
    void setForceVector(Vector3 forceVector);
private:
    Vector3 m_forceVector;
};

inline void ConstantForce::setForceVector(Vector3 forceVector) {
    m_forceVector = forceVector;
}

#endif // CONSTANTFORCE_H
