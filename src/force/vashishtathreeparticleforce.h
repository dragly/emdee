#ifndef VASHISHTATHREEPARTICLEFORCE_H
#define VASHISHTATHREEPARTICLEFORCE_H

#include <src/force/threeparticleforce.h>
#include <map>

using namespace std;

class VashishtaThreeParticleForce : public ThreeParticleForce
{
public:
    VashishtaThreeParticleForce();

    void calculateAndApplyForce(Atom *atom1, Atom *atom2, Atom *atom3);
    void calculateAndApplyForce(Atom *atom1, Atom *atom2, Atom *atom3, const Vector3 &atom2Offset, const Vector3 &atom3Offset);

    void setMapForAllPermutations(map<vector<int>, double> &theMap, const vector<int> values, double value);
    void setMapForAllPermutationsStep2(map<vector<int>, double> &theMap, const vector<int> &v, const int start, const int n, double value);
    void setMapForAllPermutationsStep2(map<vector<int>, int> &theMap, const vector<int> &v, const int start, const int n, double value);
    void setMapForAllPermutations(map<vector<int>, int> &theMap, const vector<int> values, double value);
protected:
    map<vector<int>,double> m_B;
    map<vector<int>,double> m_cosThetaBar;
    map<vector<int>,double> m_r0;
    map<vector<int>,int> m_centerAtom;
    vector<int> m_combo;

//    Vector3 force;
    Vector3 rij;
    Vector3 rik;

    Atom* atoms[3];
    Vector3 atomPosition[3];

    double force(double Bijk, double drij, double dfdrij, double drik, double dfdrik, double p, double f, double dtheta, double dpdtheta);
//    map<vector<int>,double> l;
};

#endif // VASHISHTATHREEPARTICLEFORCE_H
