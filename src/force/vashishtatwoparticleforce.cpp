#include "vashishtatwoparticleforce.h"

#include <src/atom.h>

VashishtaTwoParticleForce::VashishtaTwoParticleForce() :
    TwoParticleForce()
{
    pair<int,int> Si_Si(14,14);
    pair<int,int> Si_O(14,8);
    pair<int,int> O_Si(Si_O.second, Si_O.first);
    pair<int,int> O_O(8,8);

    H[Si_Si] = 0.057;
    H[Si_O] = 11.387;
    H[O_Si] = H[Si_O];
    H[O_O] = 51.692;
    eta[Si_Si] = 11;
    eta[Si_O] = 9;
    eta[O_Si] = eta[Si_O];
    eta[O_O] = 7;
}

void VashishtaTwoParticleForce::calculateAndApplyForce(Atom *atom1, Atom *atom2)
{
    calculateAndApplyForce(atom1, atom2, m_zeroVector);
}


void VashishtaTwoParticleForce::calculateAndApplyForce(Atom *atom1, Atom *atom2, const Vector3 &atomOffset)
{
    Vector3 rVec = atom2->position() + atomOffset - atom1->position();
    pair<int,int> combo(atom1->type().id(), atom2->type().id());
    //    cout << "Looking for " << combo.first << " " << combo.second << endl;
    double Hij = H[combo];
    double etaij = eta[combo];
    double Zi = atom1->type().effectiveCharge();
    double Zj = atom2->type().effectiveCharge();
    double alphai = atom1->type().electronicPolarizability();
    double alphaj = atom2->type().electronicPolarizability();
    double r4s = 4.43;
    double r2 = dot(rVec, rVec);
    double r = sqrt(r2);
    double r4 = r2 * r2;
//    double r5 = r2 * r2 * r;
//    double r6 = r5 * r;

    Vector3 force = -rVec*(-Hij*etaij*pow(r, -etaij)/pow(r, 2) - Zi*Zj/pow(r, 3) + (0.5*pow(Zi, 2)*alphaj + 0.5*pow(Zj, 2)*alphai)*exp(-r/r4s)/(pow(r, 5)*r4s) + 4*(0.5*pow(Zi, 2)*alphaj + 0.5*pow(Zj, 2)*alphai)*exp(-r/r4s)/pow(r, 6));

    double potential = Hij / pow(r, etaij)
            + (Zi * Zj) / r
            - (0.5 * (alphai * Zj * Zj + alphaj * Zi * Zi) / r4) * exp(-r / r4s);

    if(m_isNewtonsThirdLawEnabled) {
        atom2->addForce(force);
    }
    atom1->addForce(-force);

    if(m_isNewtonsThirdLawEnabled) {
        atom2->addPotential(0.5 * potential);
    }
    atom1->addPotential(0.5 * potential);
}
