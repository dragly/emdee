#ifndef FANNFORCE_H
#define FANNFORCE_H

#include <force/threeparticleforce.h>
struct fann;

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
    fann *m_ann;
    bool m_hasWarnedAboutMissingNetwork;
};

#endif // FANNFORCE_H
