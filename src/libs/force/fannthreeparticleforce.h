#ifndef FANNFORCE_H
#define FANNFORCE_H

#include <force/threeparticleforce.h>
struct fann;
typedef double fann_type;

class FannThreeParticleForce : public ThreeParticleForce
{
public:
    FannThreeParticleForce();

    void loadNetwork(const std::string &fileName, const std::string &boundsFilename);

    // ThreeParticleForce interface
public:
    void calculateAndApplyForce(Atom *atom1, Atom *atom2, Atom *atom3);
    void calculateAndApplyForce(Atom *atom1, Atom *atom2, Atom *atom3, const Vector3 &atom2Offset, const Vector3 &atom3Offset);

    double cutoffRadius() const;
    void setCutoffRadius(double cutoffRadius);

private:
    void warnAboutMissingNetwork();

    double m_cutoffRadius;

    fann *m_ann;
    fann_type m_fanntmp;
    bool m_hasWarnedAboutMissingNetwork;
    fann_type *testForce(fann_type *input);
};

#endif // FANNFORCE_H
