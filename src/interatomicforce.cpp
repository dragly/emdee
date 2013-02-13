#include <src/interatomicforce.h>
#include <src/moleculesystem.h>
#include <src/atom.h>

InteratomicForce::InteratomicForce() :
    tmpForce(zeros<rowvec>(3)),
    zeroVector(zeros<rowvec>(3)),
    m_potentialConstant(1),
    m_potentialConstantSquared(1),
    m_energyConstant(1),
    m_energyConstant4(4),
    m_energyConstant24(24)
{
}

void InteratomicForce::calculateAndApplyForce(Atom *atom1, Atom *atom2) {
    return calculateAndApplyForce(atom1, atom2, zeroVector);
}

void InteratomicForce::calculateAndApplyForce(Atom *atom1, Atom *atom2, const rowvec& atom2Offset)
{
    double rSquared;
    double sigmaSquaredOverRSquared;
    double sigmaSquared = m_potentialConstantSquared;
    double eps4 = m_energyConstant4;
    double eps24 = m_energyConstant24;
    tmpForce = atom2->position() + atom2Offset - atom1->position();
    // Check distances to the nearby cells
    rSquared = dot(tmpForce, tmpForce);
    sigmaSquaredOverRSquared = sigmaSquared/rSquared;
    double sigmaOverR6 = sigmaSquaredOverRSquared * sigmaSquaredOverRSquared * sigmaSquaredOverRSquared;
    double sigmaOverR12 = sigmaOverR6 * sigmaOverR6;
    // TODO Verify force term
    double factor = - ((eps24) / (rSquared)) * (2 * sigmaOverR12 - sigmaOverR6);
    tmpPotential = eps4 * (sigmaOverR12 - sigmaOverR6);

//    tmpForce *= factor; // * rVec;

    atom2->addForce(-tmpForce * factor);
    atom1->addForce(tmpForce * factor);
    atom2->addPotential(0.5 * tmpPotential);
    atom1->addPotential(0.5 * tmpPotential);
}

//const rowvec &InteratomicForce::force()
//{
//    return tmpForce;
//}

void InteratomicForce::setPotentialConstant(double potentialConstant)
{
    m_potentialConstant = potentialConstant;
    m_potentialConstantSquared = potentialConstant*potentialConstant;
}

void InteratomicForce::setEnergyConstant(double energyConstant)
{
    m_energyConstant = energyConstant;
    m_energyConstant4 = energyConstant * 4;
    m_energyConstant24 = energyConstant * 24;
}


//double InteratomicForce::potential()
//{
//    return tmpPotential;
//}
