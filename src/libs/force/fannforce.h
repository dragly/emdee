#ifndef FANNFORCE_H
#define FANNFORCE_H

#include <force/threeparticleforce.h>
struct fann;
typedef double fann_type;

class FannForce : public ThreeParticleForce
{
public:
    FannForce();

    void loadNetwork(std::string fileName);

    // ThreeParticleForce interface
public:
    void calculateAndApplyForce(Atom *atom1, Atom *atom2, Atom *atom3);
    void calculateAndApplyForce(Atom *atom1, Atom *atom2, Atom *atom3, const Vector3 &atom2Offset, const Vector3 &atom3Offset);

private:
    void warnAboutMissingNetwork();

    fann *m_ann;
    fann_type m_fanntmp;
    bool m_hasWarnedAboutMissingNetwork;
    fann_type *testForce(fann_type *input);
};

#endif // FANNFORCE_H
