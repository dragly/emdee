#include "vashishtathreeparticleforce.h"

#include <iostream>
#include <src/atom.h>
#include <stdio.h>

using namespace std;

VashishtaThreeParticleForce::VashishtaThreeParticleForce()
{
    vector<int> SiOSi({8,14,14});
    vector<int> OSiO({8,8,14});

    m_combo = vector<int>({0,0,0});
    //    B[SiOSi] = 1.4;
    //    B[OSiO] = 0.35;
    setMapForAllPermutations(m_B, SiOSi, 1.4);
    setMapForAllPermutations(m_B, OSiO, 0.35);
    setMapForAllPermutations(m_cosThetaBar, SiOSi, cos(141.0 / 360. * 2*M_PI));
    setMapForAllPermutations(m_cosThetaBar, OSiO, cos(109.47 / 360. * 2*M_PI));
    setMapForAllPermutations(m_r0, SiOSi, 2.6);
    setMapForAllPermutations(m_r0, OSiO, 2.6);
    setMapForAllPermutations(m_centerAtom, SiOSi, 8);
    setMapForAllPermutations(m_centerAtom, OSiO, 14);

    //    for(auto item : m_thetaBar) {
    //        for(int key : item.first) {
    //            cout << key << " ";
    //        }
    //        cout << item.second  << endl;
    //    }
}

/*
 * From example found on http://www.bearcave.com/random_hacks/permute.html
 */
void VashishtaThreeParticleForce::setMapForAllPermutationsStep2(map<vector<int>, double> &theMap, const vector<int> &v, const int start, const int n, double value)
{
    vector<int> v2 = v;
    if (start == n-1) {
        theMap[v2] = value;
    } else {
        for (int i = start; i < n; i++) {
            int tmp = v2[i];

            v2[i] = v2[start];
            v2[start] = tmp;
            setMapForAllPermutationsStep2(theMap, v2, start+1, n, value);
            v2[start] = v2[i];
            v2[i] = tmp;
        }
    }
}

void VashishtaThreeParticleForce::setMapForAllPermutations(map<vector<int>, double> &theMap, const vector<int> values, double value) {
    setMapForAllPermutationsStep2(theMap, values, 0, values.size(), value);
}

void VashishtaThreeParticleForce::setMapForAllPermutationsStep2(map<vector<int>, int> &theMap, const vector<int> &v, const int start, const int n, double value)
{
    vector<int> v2 = v;
    if (start == n-1) {
        theMap[v2] = value;
    } else {
        for (int i = start; i < n; i++) {
            int tmp = v2[i];

            v2[i] = v2[start];
            v2[start] = tmp;
            setMapForAllPermutationsStep2(theMap, v2, start+1, n, value);
            v2[start] = v2[i];
            v2[i] = tmp;
        }
    }
}

void VashishtaThreeParticleForce::setMapForAllPermutations(map<vector<int>, int> &theMap, const vector<int> values, double value) {
    setMapForAllPermutationsStep2(theMap, values, 0, values.size(), value);
}



void VashishtaThreeParticleForce::calculateAndApplyForce(Atom *atom1, Atom *atom2, Atom *atom3)
{
    calculateAndApplyForce(atom1, atom2, atom3, m_zeroVector, m_zeroVector);
}

void VashishtaThreeParticleForce::calculateAndApplyForce(Atom *atom1, Atom *atom2, Atom *atom3, const Vector3 &atom2Offset, const Vector3 &atom3Offset)
{
    m_combo[0] = atom1->type().number();
    m_combo[1] = atom2->type().number();
    m_combo[2] = atom3->type().number();

    if(m_centerAtom.find(m_combo) == m_centerAtom.end()) {
        return;
    }

    int centralAtom = m_centerAtom[m_combo];

    if(atom1->type().number() == centralAtom) {
        atoms[0] = atom1;
        atomPosition[0] = atom1->position();
        atoms[1] = atom2;
        atomPosition[1] = atom2->position() + atom2Offset;
        atoms[2] = atom3;
        atomPosition[2] = atom3->position() + atom3Offset;
    } else if(atom2->type().number() == centralAtom) {
        atoms[0] = atom2;
        atomPosition[0] = atom2->position() + atom2Offset;
        atoms[1] = atom1;
        atomPosition[1] = atom1->position();
        atoms[2] = atom3;
        atomPosition[2] = atom3->position() + atom3Offset;
    } else if(atom3->type().number() == centralAtom) {
        atoms[0] = atom3;
        atomPosition[0] = atom3->position() + atom3Offset;
        atoms[1] = atom1;
        atomPosition[1] = atom1->position();
        atoms[2] = atom2;
        atomPosition[2] = atom2->position() + atom2Offset;
    }

    rij = atomPosition[1] - atomPosition[0];
    rik = atomPosition[2] - atomPosition[0];

    double dotrijrik = dot(rij,rik);
    double lij2 = dot(rij, rij);
    double lik2 = dot(rik, rik);
    double lij = sqrt(lij2);
    double lik = sqrt(lik2);
    double r0 = m_r0[m_combo];
    if(lij > r0 || lik > r0) {
        return;
    }
    double shield = 1e-12; // just to avoid nan from 0 / 0 in dtheta and acos
    double l = 1.0;
    double Bijk = m_B[m_combo];
    double cosThetaBar = m_cosThetaBar[m_combo];
    double lijMinusR0 = lij - r0;
    double likMinusR0 = lik - r0;
    double lOverLijMinusR0 = l / lijMinusR0;
    double lOverLikMinusR0 = l / likMinusR0;
    double f = exp(lOverLijMinusR0 + lOverLikMinusR0);
    double dfdrij = -l*exp(lOverLikMinusR0 + lOverLijMinusR0)/(lijMinusR0*lijMinusR0);
    double dfdrik = -l*exp(lOverLikMinusR0 + lOverLijMinusR0)/(likMinusR0*likMinusR0);
    double costheta = dotrijrik / (lij * lik + shield);
    double p1 = (costheta - cosThetaBar); // TODO: Rewrite without cos and acos ;)
    double p = p1*p1;

    double sintheta = sqrt(1 - costheta*costheta);
    double dpdtheta = -2*(-cosThetaBar + costheta) * sintheta;


    double dtheta = 0;
    double drij = 0;
    double drik = 0;

    double invLij = 1 / lij;
    double invLik = 1 / lik;
    double invSqrtDotOverLenghtSquared = 1. /
                (
                sqrt(
                    1 - (dotrijrik * dotrijrik * (invLij * invLij * invLik * invLik)) + shield
                    )
                );

    if(atom1 == atoms[0] || m_isNewtonsThirdLawEnabled) {
        dtheta = 0;
        drij = 0;
        drik = 0;
        for(int a = 0; a < 3; a++) {
            dtheta = -(
                        ((dotrijrik * rik[a]) * (invLij * invLik * invLik * invLik))
                        + ((dotrijrik * rij[a]) * (invLij * invLij * invLij * invLik))
                        + ((-rij[a] - rik[a]) * (invLij * invLik))
                        ) * invSqrtDotOverLenghtSquared;
            drij = -rij[a] * invLij;
            drik = -rik[a] * invLik;
            double forceComp = force(Bijk, drij, dfdrij, drik, dfdrik, p, f, dtheta, dpdtheta);
            atoms[0]->addForce(a, forceComp);
        }
    }

    if(atom1 == atoms[1] || m_isNewtonsThirdLawEnabled) {
        dtheta = 0;
        drij = 0;
        drik = 0;
        for(int a = 0; a < 3; a++) {
            dtheta = -(
                        ((dotrijrik * (-rij[a])) * (invLij * invLij * invLij * invLik))
                        + ((rik[a]) * (invLij * invLik))
                        ) * invSqrtDotOverLenghtSquared;
            drij = rij[a] * invLij;
            double forceComp = force(Bijk, drij, dfdrij, drik, dfdrik, p, f, dtheta, dpdtheta);
            atoms[1]->addForce(a, forceComp);

        }
    }

    if(atom1 == atoms[2] || m_isNewtonsThirdLawEnabled) {
        dtheta = 0;
        drij = 0;
        drik = 0;
        for(int a = 0; a < 3; a++) {
            dtheta = -(
                        ((dotrijrik * (-rik[a])) * (invLij * invLik * invLik * invLik))
                        + ((rij[a]) * (invLij * invLik))
                        ) * invSqrtDotOverLenghtSquared;
            drik = rik[a] * invLik;
            double forceComp = force(Bijk, drij, dfdrij, drik, dfdrik, p, f, dtheta, dpdtheta);
            atoms[2]->addForce(a,forceComp);
        }
    }

    double potential = Bijk * f * p;
    atom1->addPotential(potential / 3);
    if(m_isNewtonsThirdLawEnabled) {
        atom2->addPotential(potential / 3);
        atom3->addPotential(potential / 3);
    }
}

double VashishtaThreeParticleForce::force(double Bijk, double drij, double dfdrij, double drik, double dfdrik, double p, double f, double dtheta, double dpdtheta)
{
    return -Bijk * ((drij * dfdrij + drik * dfdrik) * p + f * (dtheta * dpdtheta));
}
