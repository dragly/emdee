#include <src/interatomicforce.h>
#include <src/moleculesystem.h>
#include <src/atom.h>

InteratomicForce::InteratomicForce() :
    tmpForce(zeros<rowvec>(3)),
    zeroVector(zeros<rowvec>(3)),
    m_potentialConstant(1.2)
{
}

const rowvec InteratomicForce::force(Atom *atom1, Atom *atom2) {
    return force(atom1, atom2, zeroVector);
}

const rowvec InteratomicForce::force(Atom *atom1, Atom *atom2, const rowvec& atom2Offset)
{
    double rSquared;
    double sigmaSquaredOverRSquared;
    double sigma = m_potentialConstant;
    double eps = 1;
    rVec = atom2->position() + atom2Offset - atom1->position();
    // Check distances to the nearby cells
    rSquared = dot(rVec, rVec);
    sigmaSquaredOverRSquared = sigma*sigma/rSquared;
    double sigmaOverR6 = sigmaSquaredOverRSquared * sigmaSquaredOverRSquared * sigmaSquaredOverRSquared;
    double sigmaOverR12 = sigmaOverR6 * sigmaOverR6;
    // TODO Verify force term
    double factor = - ((24 * eps) / (rSquared)) * (2 * sigmaOverR12 - sigmaOverR6);

    tmpForce = factor * rVec;
    return tmpForce;
}

void InteratomicForce::setPotentialConstant(double potentialConstant)
{
    m_potentialConstant = potentialConstant;
}
