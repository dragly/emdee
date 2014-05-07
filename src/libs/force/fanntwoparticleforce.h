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
    double tailCorrectionMin;
    double tailCorrectionMax;
    double headCorrectionMax;
    double energyMin;
    double energyMax;
    AtomType atomType1;
    AtomType atomType2;
    double headCorrectionMaxDerivative;
    double tailCorrectionMinDerivative;
    double headCorrectionMaxEnergy;
    double tailCorrectionMinEnergy;
    double energyOffset;

    double rescaleDistance(double r12) const;

    double rescaleEnergy(double energy) const;

    double rescaleEnergyDerivative(double value) const;

    double tailCorrectionEnergy(double l12) const;
    double tailCorrectionDerivative(double l12) const;
    double headCorrectionEnergy(double l12) const;
    double headCorrectionForce(double l12) const;
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
