#include "vashishtathreeparticleforce.h"

#include <iostream>
#include <src/atom.h>

using namespace std;

VashishtaThreeParticleForce::VashishtaThreeParticleForce()
{
}

void VashishtaThreeParticleForce::calculateAndApplyForce(Atom *atom1, Atom *atom2, Atom *atom3)
{
    calculateAndApplyForce(atom1, atom2, atom3, m_zeroVector, m_zeroVector);
}

void VashishtaThreeParticleForce::calculateAndApplyForce(Atom *atomi, Atom *atomj, Atom *atomk, const Vector3 &atomjOffset, const Vector3 &atomkOffset)
{
    Vector3 rij = atomj->position() + atomjOffset - atomi->position();
    Vector3 rik = atomk->position() + atomkOffset - atomk->position();
    double dotrijrik = dot(rij,rik);
    double lij2 = dot(rij, rij);
    double lik2 = dot(rik, rik);
    double lij = sqrt(lij2);
    double lik = sqrt(lik2);
    double l = 1.0;
    double Bijk = 1.0;
    double thetabar = 109;
    double r0 = 2.6;

    double theta = acos(dotrijrik / (lij * lik));

    double f;
    if(lij < r0 || lik < r0) {
        f = exp(l / (lij - r0) + l / (lik - r0));
    } else {
        f = 0;
    }
    double p1 = (cos(theta) - cos(thetabar)); // TODO: Rewrite without cos and acos ;)
    double p = p1*p1;

    double dfdrij = -2*l*exp(2*l/(-r0 + lij))/pow(-r0 + lij, 2);
    double dfdrik = -l*exp(l/(-r0 + lik) + l/(-r0 + lij))/pow(-r0 + lik, 2);
//    double dfdtheta = 0;

//    double dpdrij = 0;
//    double dpdrik = 0;
    double dpdtheta = -2*(-cos(thetabar) + cos(theta)) * sin(theta);

    for(int iAtom = 0; iAtom < 3; iAtom++) {
        Vector3 force;
        for(int a = 0; a < 3; a++) {
            double dtheta = 0;

            switch(iAtom) {
            case 0:
                dtheta = (
                            ((dotrijrik * rik[a]) / (lij * lik * lik * lik))
                            + ((dotrijrik * rij[a]) / (lij * lij * lij * lik))
                            + ((-rij[a] - rik[a]) / (lij * lik))
                            ) / (
                            sqrt(
                                1 - (dotrijrik * dotrijrik / (lij * lij * lik * lik))
                                )
                            );
                break;
            case 1:
                dtheta = (
                            ((dotrijrik * (-rij[a])) / (lij * lik * lik * lik))
                            + ((rik[a]) / (lij * lik))
                            ) / (
                            sqrt(
                                1 - (dotrijrik * dotrijrik / (lij * lij * lik * lik))
                                )
                            );

                break;
            case 2:
                dtheta = (
                            ((dotrijrik * (-rik[a])) / (lij * lik * lik * lik))
                            + ((rij[a]) / (lij * lik))
                            ) / (
                            sqrt(
                                1 - (dotrijrik * dotrijrik / (lij * lij * lik * lik))
                                )
                            );
                break;
            }

            double drij = 0;
            if(iAtom == 0) {
                drij = -rij[a] / lij;
            } else if(iAtom == 1) {
                drij = rij[a] / lij;
            }

            double drik = 0;
            if(iAtom == 0) {
                drik = -rik[a] / lik;
            } else if(iAtom == 1) {
                drik = rik[a] / lik;
            }

            force[a] = Bijk * ((drij * dfdrij + drik * dfdrik) * p + f * (dtheta * dpdtheta));
        }
        switch(iAtom) {
        case 0:
            atomi->addForce(force);
            break;
        case 1:
            atomj->addForce(force);
            break;
        case 2:
            atomk->addForce(force);
            break;
        }
    }
}
