#ifndef FANNFORCE_H
#define FANNFORCE_H

#include <force/threeparticleforce.h>
struct fann;
typedef double fann_type;

class FannThreeParticleNetwork {
public:
    fann* ann;
    double r12Min;
    double r12Max;
    double r13Min;
    double r13Max;
    double minDistance;
    double energyMin;
    double energyMax;
    AtomType atomType1;
    AtomType atomType2;
    AtomType atomType3;
    double headCorrectionMaxDerivative;
    double tailCorrectionMinDerivative;
    double headCorrectionMaxEnergy;
    double tailCorrectionMinEnergy;
    double energyOffset;
    double angleMin = M_PI / 10;
    double angleMax = M_PI;

    double rescaleDistance12(double r12) const;
    double rescaleDistance13(double r12) const;
    double rescaleEnergy(double energy) const;
//    double rescaleEnergyDerivative(double value) const;
//    double tailCorrectionDampingFactor(double l12) const;
//    double tailCorrectionDampingFactorDerivative(double l12) const;
//    double headCorrectionEnergy(double l12) const;
//    double headCorrectionForce(double l12) const;
    double rescaleEnergyDerivativeR12(double value) const;
    double rescaleEnergyDerivativeR13(double value) const;
    double rescaleEnergyDerivativeAngle(double value) const;
};

class FannThreeParticleForce : public ThreeParticleForce
{
public:
    FannThreeParticleForce();

    void loadNetwork(const std::string &fileName, const std::string &boundsFilename, double minDistance = 0.0);

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
    vector<FannThreeParticleNetwork> m_networks;
};

#endif // FANNFORCE_H
