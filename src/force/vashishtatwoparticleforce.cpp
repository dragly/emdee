#include "vashishtatwoparticleforce.h"

#include <src/atom.h>

VashishtaTwoParticleForce::VashishtaTwoParticleForce() :
    nParticleTypes(2),
    TwoParticleForce(),
    invsqrtpi(1 / sqrt(M_PI))
{
//    pair<int,int> Si_Si(14,14);
//    pair<int,int> Si_O(14,8);
//    pair<int,int> O_Si(Si_O.second, Si_O.first);
//    pair<int,int> O_O(8,8);
    int nPermutations = nParticleTypes*nParticleTypes;
    H = zeros(nPermutations);
    eta = zeros(nPermutations);

    int Si_Si[2] = {0,0};
    int Si_O[2] = {0,1};
    int O_Si[2] = {1,0};
    int O_O[2] = {1,1};

    H[comboHash(Si_Si)] = 0.057;
    H[comboHash(Si_O)] = 11.387;
    H[comboHash(O_Si)] = H[comboHash(Si_O)];
    H[comboHash(O_O)] = 51.692;
    eta[comboHash(Si_Si)] = 11;
    eta[comboHash(Si_O)] = 9;
    eta[comboHash(O_Si)] = eta[comboHash(Si_O)];
    eta[comboHash(O_O)] = 7;
}

int VashishtaTwoParticleForce::comboHash(int v[2]) {
    return v[0] + v[1] * nParticleTypes;
}

void VashishtaTwoParticleForce::calculateAndApplyForce(Atom *atom1, Atom *atom2)
{
    calculateAndApplyForce(atom1, atom2, m_zeroVector);
}


void VashishtaTwoParticleForce::calculateAndApplyForce(Atom *atom1, Atom *atom2, const Vector3 &atomOffset)
{
    rVec = atom2->position() + atomOffset - atom1->position();
//    pair<int,int> combo(atom1->type().number(), atom2->type().number());
    int v[2];
    v[0] = atom1->type().index();
    v[1] = atom2->type().index();
    int combo = comboHash(v);
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

    double factorH = -Hij*etaij*pow(r, -etaij) * invr2;
    double factorExp = (0.5*Zi2*alphaj + 0.5*Zj2*alphai)*expR4OverR4s*invr5*invr4s + 4*(0.5*Zi2*alphaj + 0.5*Zj2*alphai)*expR4OverR4s * invr6;

    double ewaldAlpha = 1 / r4s; // From Toukmaji, Board, Computer Physics Communications 95 (1996)
    double factorCoulomb = -Zi * Zj * invr3 * (erfc(ewaldAlpha*r) + 2*ewaldAlpha * invsqrtpi * r * exp(-(ewaldAlpha * r)*(ewaldAlpha * r)));
//    double factorCoulomb =  - Zi*Zj * invr3;

    force = -rVec*(factorH + factorExp + factorCoulomb);

    double potentialH = Hij * pow(invr, etaij);
    double potentialCoulomb = (Zi * Zj) * invr;
    double potentialExp = - (0.5 * (alphai * Zj * Zj + alphaj * Zi * Zi) * invr4) * expR4OverR4s;

    double potential = potentialH + potentialCoulomb + potentialExp;

    if(m_isNewtonsThirdLawEnabled) {
        atom2->addForce(force);
    }
    atom1->addForce(-force);

    if(m_isNewtonsThirdLawEnabled) {
        atom2->addPotential(0.5 * potential);
    }
    atom1->addPotential(0.5 * potential);
}
