#include "vashishtatwoparticleforce.h"

#include <src/atom.h>

VashishtaTwoParticleForce::VashishtaTwoParticleForce() :
    TwoParticleForce()
{
}

void VashishtaTwoParticleForce::calculateAndApplyForce(Atom *atom1, Atom *atom2)
{
    calculateAndApplyForce(atom1, atom2, m_zeroVector);
}


void VashishtaTwoParticleForce::calculateAndApplyForce(Atom *atom1, Atom *atom2, const Vector3 &atomOffset)
{
    Vector3 rVec = atom2->position() + atomOffset - atom1->position();
    double Hij = 1.7;
    double etaij = 1.6;
    double Z_i = 1.5;
    double Z_j = 1.4;
    double alphai = 1.3;
    double alphaj = 1.2;
    double r4s = 1.1;
    double r2 = dot(rVec, rVec);
    double r = sqrt(r2);
    double r5 = r2 * r2 * r;
    double r6 = r5 * r;


    Vector3 force =rVec*(-1.0*Hij*etaij*pow(r, -etaij)/pow(r, 2) - 0.5*pow(Z_i, 2)*alphaj*exp(-r/r4s)/(r5*r4s) - 2.0*pow(Z_i, 2)*alphaj*exp(-r/r4s)/r6 - 1.0*Z_i*Z_j/pow(r, 3) - 0.5*pow(Z_j, 2)*alphai*exp(-r/r4s)/(r5*r4s) - 2.0*pow(Z_j, 2)*alphai*exp(-r/                                                                                              r4s)/r6);

    if(m_isNewtonsThirdLawEnabled) {
        atom2->addForce(force);
    }
    atom1->addForce(-force);

}
