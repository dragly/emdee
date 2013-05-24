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
    pair<int,int> combo(atom1->type().number(), atom2->type().number());
    //    cout << "Looking for " << combo.first << " " << combo.second << endl;
    double Hij = H[combo];
    double etaij = eta[combo];
    double Zi = atom1->type().effectiveCharge();
    double Zj = atom2->type().effectiveCharge();
    double Zi2 = Zi * Zi;
    double Zj2 = Zj * Zj;
    double alphai = atom1->type().electronicPolarizability();
    double alphaj = atom2->type().electronicPolarizability();
    double r4s = 4.43;
    double r2 = dot(rVec, rVec);
    double r = sqrt(r2);
//    double r3 = r2 * r;
//    double r4 = r2 * r2;
//    double r5 = r2 * r2 * r;
//    double r6 = r5 * r;
    double invr = 1 / r;
    double invr2 = invr * invr;
    double invr3 = invr2 * invr;
    double invr4 = invr2 * invr2;
    double invr5 = invr4 * invr;
    double invr6 = invr3 * invr3;
    double invr4s = 1 / r4s;
    double expR4OverR4s = exp(-r/r4s);

    Vector3 force = -rVec*(-Hij*etaij*pow(r, -etaij) * invr2 - Zi*Zj * invr3 + (0.5*Zi2*alphaj + 0.5*Zj2*alphai)*expR4OverR4s*invr5*invr4s + 4*(0.5*Zi2*alphaj + 0.5*Zj2*alphai)*expR4OverR4s * invr6);

    double potential = Hij * pow(invr, etaij)
            + (Zi * Zj) * invr
            - (0.5 * (alphai * Zj * Zj + alphaj * Zi * Zi) * invr4) * expR4OverR4s;

    if(m_isNewtonsThirdLawEnabled) {
        atom2->addForce(force);
    }
    atom1->addForce(-force);

    if(m_isNewtonsThirdLawEnabled) {
        atom2->addPotential(0.5 * potential);
    }
    atom1->addPotential(0.5 * potential);
}
