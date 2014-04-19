#ifndef FANNTWOPARTICLEFORCE_H
#define FANNTWOPARTICLEFORCE_H

#include <vector>
#include "atomtype.h"
#include <force/twoparticleforce.h>
struct fann;
typedef double fann_type;

using std::vector;

class FannTwoParticleNetwork {
public:
    fann* ann;
    double r12Min;
    double r12Max;
    double energyMin;
    double energyMax;
    AtomType atomType1;
    AtomType atomType2;

    double rescaleDistance(double r12) const
    {
        return (r12 - r12Min) / (r12Max - r12Min) * 0.8 + 0.1;
    }

    double rescaleEnergy(double energy) const
    {
        return ((energy - 0.1) / 0.8) * (energyMax - energyMin) + energyMin;
    }
};

class FannTwoParticleForce : public TwoParticleForce
{
public:
    FannTwoParticleForce();

    void addNetwork(const AtomType& atomType1, const AtomType& atomType2, const std::string& fileName, const std::string& boundsFilename);

    // TwoParticleForce interface
public:
    virtual void calculateAndApplyForce(Atom *atom1, Atom *atom2);
    virtual void calculateAndApplyForce(Atom *atom1, Atom *atom2, const Vector3 &atom2Offset);

private:
    void warnAboutMissingNetwork();
    fann_type m_fanntmp;
    bool m_hasWarnedAboutMissingNetwork;
    vector<FannTwoParticleNetwork> m_networks;
};

#endif // FANNTWOPARTICLEFORCE_H
