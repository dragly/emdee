#ifndef FANNTWOPARTICLEFORCE_H
#define FANNTWOPARTICLEFORCE_H

#include <force/twoparticleforce.h>
struct fann;
typedef double fann_type;

class FannTwoParticleForce : public TwoParticleForce
{
public:
    FannTwoParticleForce();

    void loadNetwork(const std::string& fileName, const std::string& boundsFilename);

    // TwoParticleForce interface
public:
    virtual void calculateAndApplyForce(Atom *atom1, Atom *atom2);
    virtual void calculateAndApplyForce(Atom *atom1, Atom *atom2, const Vector3 &atom2Offset);

private:
    void warnAboutMissingNetwork();
    double rescaleDistance(double r12);
    double rescaleEnergy(double energy);
    fann *m_ann;
    fann_type m_fanntmp;
    bool m_hasWarnedAboutMissingNetwork;
    fann_type *testForce(fann_type *input);

    double m_r12Min;
    double m_r12Max;
    double m_energyMin;
    double m_energyMax;
};

#endif // FANNTWOPARTICLEFORCE_H
