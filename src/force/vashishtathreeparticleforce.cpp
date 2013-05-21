#include "vashishtathreeparticleforce.h"

#include <iostream>
#include <src/atom.h>
#include <stdio.h>

using namespace std;

VashishtaThreeParticleForce::VashishtaThreeParticleForce()
{
    vector<int> SiOSi({8,14,14});
    vector<int> OSiO({8,8,14});
//    B[SiOSi] = 1.4;
//    B[OSiO] = 0.35;
    setMapForAllPermutations(m_B, SiOSi, 1.4);
    setMapForAllPermutations(m_B, OSiO, 0.35);
    setMapForAllPermutations(m_thetaBar, SiOSi, 141.0 / 360. * 2*M_PI);
    setMapForAllPermutations(m_thetaBar, OSiO, 109.47 / 360. * 2*M_PI);
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
    vector<int> combo({atom1->type().id(), atom2->type().id(), atom3->type().id()});
    if(m_centerAtom.find(combo) == m_centerAtom.end()) {
        return;
    }

    int centralAtom = m_centerAtom[combo];

    Atom* atomi;
    Atom* atomj;
    Atom* atomk;

    Vector3 atomiPosition;
    Vector3 atomjPosition;
    Vector3 atomkPosition;

    if(atom1->type().id() == centralAtom) {
        atomi = atom1;
        atomiPosition = atom1->position();
        atomj = atom2;
        atomjPosition = atom2->position() + atom2Offset;
        atomk = atom3;
        atomkPosition = atom3->position() + atom3Offset;
    } else if(atom2->type().id() == centralAtom) {
        atomi = atom2;
        atomiPosition = atom2->position() + atom2Offset;
        atomj = atom1;
        atomjPosition = atom1->position();
        atomk = atom3;
        atomkPosition = atom3->position() + atom3Offset;
    } else if(atom3->type().id() == centralAtom) {
        atomi = atom3;
        atomiPosition = atom3->position() + atom3Offset;
        atomj = atom1;
        atomjPosition = atom1->position();
        atomk = atom2;
        atomkPosition = atom2->position() + atom2Offset;
    }

    Vector3 rij = atomjPosition - atomiPosition;
    Vector3 rik = atomkPosition - atomiPosition;
    double dotrijrik = dot(rij,rik);
    double lij2 = dot(rij, rij);
    double lik2 = dot(rik, rik);
    double lij = sqrt(lij2);
    double lik = sqrt(lik2);
    double l = 1.0;
    double Bijk = m_B[combo];
    double thetabar = m_thetaBar[combo];
    double r0 = m_r0[combo];

    double theta = acos(dotrijrik / (lij * lik));

    //    cout << "theta " << theta << endl;

    //    cout << rij << " " << rik << endl;
    //    cout << "dotrijrik " << dotrijrik << endl;
    //    cout << "theta " << theta << endl;
    //    cout << "acos(dotrijrik / (lij * lik)) =" << "acos(" << dotrijrik << "/" << "(" << lij << "*" << lik << "))";

    double f = 0;
    double dfdrij = 0;
    double dfdrik = 0;
    if(lij < r0 && lik < r0) {
        f = exp(l / (lij - r0) + l / (lik - r0));
        dfdrij = -l*exp(l/(-r0 + lik) + l/(-r0 + lij))/pow(-r0 + lij, 2);
        dfdrik = -l*exp(l/(-r0 + lik) + l/(-r0 + lij))/pow(-r0 + lik, 2);
    }

    double p1 = (cos(theta) - cos(thetabar)); // TODO: Rewrite without cos and acos ;)
    double p = p1*p1;

    //    cout << "p " << p << endl;
    //    cout << "f " << f << endl;


    //    double dfdtheta = 0;

    //    double dpdrij = 0;
    //    double dpdrik = 0;
    double dpdtheta = -2*(-cos(thetabar) + cos(theta)) * sin(theta);

    for(int iAtom = 0; iAtom < 3; iAtom++) {
        Atom* currentAtom;
        switch(iAtom) {
        case 0:
            currentAtom = atomi;
            break;
        case 1:
            currentAtom = atomj;
            break;
        case 2:
            currentAtom = atomk;
            break;
        }

        if(!m_isNewtonsThirdLawEnabled && currentAtom != atom1) {
            continue;
        }
        Vector3 force;
        for(int a = 0; a < 3; a++) {
            double dtheta = 0;

            double shield = 1e-12; // just to avoid nan from 0 / 0 in dtheta
            switch(iAtom) {
            case 0:
                dtheta = -(
                            ((dotrijrik * rik[a]) / (lij * lik * lik * lik))
                            + ((dotrijrik * rij[a]) / (lij * lij * lij * lik))
                            + ((-rij[a] - rik[a]) / (lij * lik))
                            ) / (
                            sqrt(
                                1 - (dotrijrik * dotrijrik / (lij * lij * lik * lik)) + shield
                                )
                            );
                break;
            case 1:
                dtheta = -(
                            ((dotrijrik * (-rij[a])) / (lij * lij * lij * lik))
                            + ((rik[a]) / (lij * lik))
                            ) / (
                            sqrt(
                                1 - (dotrijrik * dotrijrik / (lij * lij * lik * lik)) + shield
                                )
                            );

                break;
            case 2:
                dtheta = -(
                            ((dotrijrik * (-rik[a])) / (lij * lik * lik * lik))
                            + ((rij[a]) / (lij * lik))
                            ) / (
                            sqrt(
                                1 - (dotrijrik * dotrijrik / (lij * lij * lik * lik)) + shield
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
            } else if(iAtom == 2) {
                drik = rik[a] / lik;
            }

            force[a] = -Bijk * ((drij * dfdrij + drik * dfdrik) * p + f * (dtheta * dpdtheta));
            a;
            //            force[a] = Bijk * (drij * dfdrij + drik * dfdrik);
            //            force[a] = Bijk * dtheta * dpdtheta;
        }

        currentAtom->addForce(force);
    }

    double potential = Bijk * f * p;
    //    double potential = Bijk * f;
    //    double potential = Bijk * p;
    atomi->addPotential(potential / 3);
    if(m_isNewtonsThirdLawEnabled) {
        atomj->addPotential(potential / 3);
        atomk->addPotential(potential / 3);
    }
}
