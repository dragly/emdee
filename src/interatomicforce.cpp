#include "interatomicforce.h"
#include "moleculesystem.h"
#include "atom.h"

InteratomicForce::InteratomicForce(MoleculeSystem *parent) :
    tmpForce(zeros<rowvec>(3)),
    zeroVector(zeros<rowvec>(3)),
    moleculeSystem(parent)
{
}

const rowvec InteratomicForce::force(Atom *atom1, Atom *atom2) {
    return force(atom1, atom2, zeroVector);
}

const rowvec InteratomicForce::force(Atom *atom1, Atom *atom2, const rowvec& atom2Offset)
{
    double rSquared;
    double sigmaSquaredOverRSquared;
    double sigma = moleculeSystem->potentialConstant();
    double eps = 1;
    rVec = atom2->position() + atom2Offset - atom1->position();
    // Check distances to the nearby cells
    rSquared = dot(rVec, rVec);
    sigmaSquaredOverRSquared = sigma*sigma/rSquared;
    // TODO Verify force term
    double factor = - ((24 * eps) / (rSquared)) * (2 * pow((sigmaSquaredOverRSquared), 6) - pow((sigmaSquaredOverRSquared), 3));

    tmpForce = factor * rVec;
    return tmpForce;
}
