#ifndef VASHISHTATWOPARTICLEFORCE_H
#define VASHISHTATWOPARTICLEFORCE_H

#include <src/force/twoparticleforce.h>

#include <map>

using namespace std;

class VashishtaTwoParticleForce : public TwoParticleForce
{
public:
    VashishtaTwoParticleForce();

    void calculateAndApplyForce(Atom *atom1, Atom *atom2);
    void calculateAndApplyForce(Atom *atom1, Atom *atom2, const Vector3 &atomOffset);

//    map<pair<int,int>, double> H;
//    map<pair<int,int>, double> eta;

    vec H;
    vec eta;

    int nParticleTypes;
    int comboHash(int v[]);
    Vector3 rVec;
    Vector3 force;
    double invsqrtpi;
};

#endif // VASHISHTATWOPARTICLEFORCE_H
