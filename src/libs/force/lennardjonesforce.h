#ifndef LENNARDJONESFORCE_H
#define LENNARDJONESFORCE_H

#include <force/twoparticleforce.h>

class LennardJonesForce : public TwoParticleForce
{
public:
    LennardJonesForce();

    void calculateAndApplyForce(Atom *atom1, Atom *atom2);
    void calculateAndApplyForce(Atom* atom1, Atom* atom2, const Vector3 &atom2Offset);

    void setPotentialConstant(double potentialConstant);
    void setEnergyConstant(double energyConstant);
    void setShift(double shift);

protected:
    Vector3 zeroVector;
    Vector3 tmpVector;

    double m_potentialConstant;
    double m_potentialConstantSquared;
    double m_energyConstant;
    double m_energyConstant4;
    double m_energyConstant24;
    double m_shift;
    double m_shiftSquared;
};

#endif // LENNARDJONESFORCE_H
