#include "lennardjonesforce.h"
#include <src/atom.h>

LennardJonesForce::LennardJonesForce() :
    //    tmpForce(zeros<rowvec>(3)),
    zeroVector(zeros<rowvec>(3)),
    m_potentialConstant(1),
    m_potentialConstantSquared(1),
    m_energyConstant(1),
    m_energyConstant4(4),
    m_energyConstant24(24)
{
}


void LennardJonesForce::calculateAndApplyForce(Atom *atom1, Atom *atom2) {
    return calculateAndApplyForce(atom1, atom2, zeroVector);
}

void LennardJonesForce::calculateAndApplyForce(Atom *atom1, Atom *atom2, const rowvec& atom2Offset)
{
    double rSquared;
    double sigmaSquaredOverRSquared;
    double &sigmaSquared = m_potentialConstantSquared;
    double &eps4 = m_energyConstant4;
    double &eps24 = m_energyConstant24;
    // Check distances to the nearby cells
    double x = atom2->position()(0) + atom2Offset(0) - atom1->position()(0);
    double y = atom2->position()(1) + atom2Offset(1) - atom1->position()(1);
    double z = atom2->position()(2) + atom2Offset(2) - atom1->position()(2);
    rSquared = x*x + y*y + z*z;
    sigmaSquaredOverRSquared = sigmaSquared/rSquared;
    double sigmaOverR6 = sigmaSquaredOverRSquared * sigmaSquaredOverRSquared * sigmaSquaredOverRSquared;
    double sigmaOverR12 = sigmaOverR6 * sigmaOverR6;
    double factor = ((eps24) / (rSquared)) * (2 * sigmaOverR12 - sigmaOverR6);

    // The following is the same as
    //    atom2->m_force -= tmpForce * factor;
    //    atom2->m_parent->m_force -= tmpForce * factor;
    //    atom1->m_force += tmpForce * factor;
    //    atom1->m_parent->m_force += tmpForce * factor;
    //    atom1->m_force += tmpForce * factor;

    // Acting directly on elements
    atom2->addForce(0, x * factor);
    atom2->addForce(1, y * factor);
    atom2->addForce(2, z * factor);
    atom1->addForce(0, -x * factor);
    atom1->addForce(1, -y * factor);
    atom1->addForce(2, -z * factor);

    // Potential
    double tmpPotential = eps4 * (sigmaOverR12 - sigmaOverR6);
    atom2->addPotential(0.5 * tmpPotential);
    atom1->addPotential(0.5 * tmpPotential);

    // Pressure
    double pressure = factor * rSquared; // dot product
    atom2->addLocalPressure(0.5 * pressure);
    atom1->addLocalPressure(0.5 * pressure);

}

void LennardJonesForce::setPotentialConstant(double potentialConstant)
{
    m_potentialConstant = potentialConstant;
    m_potentialConstantSquared = potentialConstant*potentialConstant;
}

void LennardJonesForce::setEnergyConstant(double energyConstant)
{
    m_energyConstant = energyConstant;
    m_energyConstant4 = energyConstant * 4;
    m_energyConstant24 = energyConstant * 24;
}