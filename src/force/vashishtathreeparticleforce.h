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
protected:
    map<vector<int>,double> m_B;
    map<vector<int>,double> m_thetaBar;
    map<vector<int>,double> m_r0;
//    map<vector<int>,double> l;
};

#endif // VASHISHTATHREEPARTICLEFORCE_H
