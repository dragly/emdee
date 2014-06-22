#include "lennardjonesforce.h"
#include <atom.h>

/*!
 * \class LennardJonesForce
 * \brief The LennardJonesForce class calculates the force and potential
 * based on the Lennard-Jones potential.
 */

LennardJonesForce::LennardJonesForce() :
    m_potentialConstant(1),
    m_potentialConstantSquared(1),
    m_energyConstant(1),
    m_energyConstant4(4),
    m_energyConstant24(24),
    m_shift(0.0)
{
}


void LennardJonesForce::calculateAndApplyForce(Atom *atom1, Atom *atom2) {
    TwoParticleForce::calculateAndApplyForce(atom1, atom2);
}

void LennardJonesForce::calculateAndApplyForce(Atom *atom1, Atom *atom2, const Vector3& atom2Offset)
{
    double rSquared;
    double sigmaSquaredOverRSquared;
    double &sigmaSquared = m_potentialConstantSquared;
    double &eps4 = m_energyConstant4;
    double &eps24 = m_energyConstant24;
    // Check distances to the nearby cells

    //    tmpVector += atom2->position() + atom2Offset - atom1->position();
    //    tmpVector = atom2->position();
    //    tmpVector += atom2Offset;
    //    tmpVector -= atom1->position();
    // For some silly reason this is faster than the line above...
    double x;
    double y;
    double z;
    x = atom2->position()(0) + atom2Offset(0) - atom1->position()(0);
    y = atom2->position()(1) + atom2Offset(1) - atom1->position()(1);
    z = atom2->position()(2) + atom2Offset(2) - atom1->position()(2);
    rSquared = x*x + y*y + z*z;
    double r = sqrt(rSquared) + m_shift;
    rSquared = r*r;
    sigmaSquaredOverRSquared = sigmaSquared/rSquared;
    double sigma6OverR6 = sigmaSquaredOverRSquared * sigmaSquaredOverRSquared * sigmaSquaredOverRSquared;
    double sigma12OverR12 = sigma6OverR6 * sigma6OverR6;
    double factor = ((eps24) / (rSquared)) * (2 * sigma12OverR12 - sigma6OverR6);

    // The following is the same as
    //    atom2->m_force -= tmpForce * factor;
    //    atom2->m_parent->m_force -= tmpForce * factor;
    //    atom1->m_force += tmpForce * factor;
    //    atom1->m_parent->m_force += tmpForce * factor;
    //    atom1->m_force += tmpForce * factor;

    // Acting directly on elements
    if(isNewtonsThirdLawEnabled()) {
        atom2->addForce(0, x * factor);
        atom2->addForce(1, y * factor);
        atom2->addForce(2, z * factor);
    }
    atom1->addForce(0, -x * factor);
    atom1->addForce(1, -y * factor);
    atom1->addForce(2, -z * factor);

    if(isCalculatePotentialEnabled()) {
        // Potential
        double tmpPotential = eps4 * (sigma12OverR12 - sigma6OverR6);
        if(isNewtonsThirdLawEnabled()) {
            atom2->addPotential(0.5 * tmpPotential);
        }
        atom1->addPotential(0.5 * tmpPotential);
    }

    if(isCalculatePressureEnabled()) {
        // Virial pressure on particle 1 from 2
        double pressure = factor * rSquared; // Fij * rij
        // For better visualization, half the pressure is added to atom1 and atom2
        atom1->addLocalPressure(0.5 * pressure);
        if(isNewtonsThirdLawEnabled()) {
            atom2->addLocalPressure(0.5 * pressure);
        }
    }
}

void LennardJonesForce::setPotentialConstant(double potentialConstant)
{
    m_potentialConstant = potentialConstant;
    m_potentialConstantSquared = potentialConstant*potentialConstant;
    setCutoffRadius(3*potentialConstant);
}

void LennardJonesForce::setEnergyConstant(double energyConstant)
{
    m_energyConstant = energyConstant;
    m_energyConstant4 = energyConstant * 4;
    m_energyConstant24 = energyConstant * 24;
}

void LennardJonesForce::setShift(double shift)
{
    m_shift = shift;
    m_shiftSquared = shift*shift;
}
