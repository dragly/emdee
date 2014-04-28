#include "fannthreeparticleforce.h"

#include <atom.h>
#include <math/vector3.h>
#include <utils/logging.h>
#include <utils/fannderivative.h>

#include <doublefann.h>
#include <iomanip>

double l12Min = 1.0;
double l12Max = 5.0;
double l13Min = 1.0;
double l13Max = 5.0;
double angleMin = M_PI / 10;
double angleMax = M_PI;
double energyMin = 0.0;
double energyMax = 1.0;
double energyOffset = 0.0;

double rescale(double value, double valueMin, double valueMax) {
    return (value - valueMin) / (valueMax - valueMin) * 0.8 + 0.1;
}

double rescaleEnergy(double value, double valueMin, double valueMax) {
    return ((value - 0.1) / 0.8) * (valueMax - valueMin) + valueMin;
}

double rescaleEnergyDerivative(double derivative, double valueMin, double valueMax, double energyMin, double energyMax)
{
    return (energyMax - energyMin) / (valueMax - valueMin) * derivative;
}

FannThreeParticleForce::FannThreeParticleForce() :
    m_ann(0),
    m_hasWarnedAboutMissingNetwork(false)
{
}

void FannThreeParticleForce::warnAboutMissingNetwork()
{
    if(!m_hasWarnedAboutMissingNetwork) {
        LOG(ERROR) << "FANN network not loaded. Cannot apply force." << endl;
        m_hasWarnedAboutMissingNetwork = true;
    }
}
double FannThreeParticleForce::cutoffRadius() const
{
    return m_cutoffRadius;
}

void FannThreeParticleForce::setCutoffRadius(double cutoffRadius)
{
    m_cutoffRadius = cutoffRadius;
}


void FannThreeParticleForce::loadNetwork(const std::string& fileName,
                                         const std::string& boundsFilename)
{
    m_ann = fann_create_from_file(fileName.c_str());

    // Read bounds
    ifstream boundsFile(boundsFilename);
    int nInputs;
    boundsFile >> nInputs;
    int nOutputs;
    boundsFile >> nOutputs;

    boundsFile >> l12Min;
    boundsFile >> l12Max;
    boundsFile >> l13Min;
    boundsFile >> l13Max;
    boundsFile >> angleMin;
    boundsFile >> angleMax;
    boundsFile >> energyMin;
    boundsFile >> energyMax;

    cout << "3P R12 bounds: " << l12Min << "," << l12Max << endl;
    cout << "3P R13 bounds: " << l13Min << "," << l13Max << endl;
    cout << "3P Energy bounds: " << energyMin << "," << energyMax << endl;
}

void FannThreeParticleForce::calculateAndApplyForce(Atom *atom1, Atom *atom2, Atom *atom3)
{
    calculateAndApplyForce(atom1, atom2, atom3, m_zeroVector, m_zeroVector);
}

//fann_type* FannThreeParticleForce::testForce(fann_type* input) {
//    m_fanntmp = 0;

//    // Scaling from ångstrøm to atomic units
//    double siToAU = 1.8897;
//    double eqDistanceScaled = rescale(1.0*siToAU, l12Min, l12Max);
//    double eqAngleScaled = rescale(1.9024, angleMin, angleMax);
//    double k1 = 50.0;
//    double k2 = 1.0;
//    m_fanntmp += k1*(input[0] - eqDistanceScaled) * (input[0] - eqDistanceScaled);
//    m_fanntmp += k1*(input[1] - eqDistanceScaled) * (input[1] - eqDistanceScaled);
//    m_fanntmp += k2*(input[2] - eqAngleScaled) * (input[2] - eqAngleScaled);
//    //    cout << input[2] << " vs. " << eqAngleScaled << endl;
//    return &m_fanntmp;
//}

void FannThreeParticleForce::calculateAndApplyForce(Atom *atom1, Atom *atom2, Atom *atom3, const Vector3 &atom2Offset, const Vector3 &atom3Offset)
{
    bool symmetric = true;
    bool damping = true;
    if(!m_ann) {
        warnAboutMissingNetwork();
        return;
    }

    bool equalParticles = (atom1->type() == atom2->type() && atom1->type() == atom3->type());
    bool atom1HasLowestID = (atom1->id() < atom2->id() && atom1->id() < atom3->id());

    if(equalParticles && !atom1HasLowestID) {
        // Return if the atom types are equal and atom1 doesn't have the lowest ID
        return;
    }

    //    if(!(atom1->type().number() == 8 && atom2->type().number() == 1 && atom3->type().number() == 1)) {
    //        return;
    //    }
    // TODO Make this general
    if(!(atom1->type().number() == 1 && atom2->type().number() == 1 && atom3->type().number() == 1)) {
        return;
    }

    Vector3 r12 = atom2->position() + atom2Offset - atom1->position();
    Vector3 r13 = atom3->position() + atom3Offset - atom1->position();
    double dotr12r13 = dot(r12, r13);
    double l12Squared = dot(r12, r12);
    double l13Squared = dot(r13, r13);
    double l12 = sqrt(l12Squared);
    double l13 = sqrt(l13Squared);
    // TODO: Use cos angle as parameter instead of angle
    //    double angle2 = acos((l12*l12 + l13*l13 - l23*l23) / (2 * l12 * l13));
    double angleParam = dotr12r13 / (l12*l13);
    if(angleParam > 1.0 || angleParam < -1.0) {
        // This should never happen, but means that atom2 and atom3 are on top of each other
        return;
    }
    double angle = acos(angleParam);
    if(angle < M_PI / 3.0) {
        // If the angle is less than 60 degrees (pi/3), we should swap the atoms so that we get one with
        // higher angle. This way we don't have to train the neural network for angles below pi/3
        // Swapping with the closest atom should give the highest angle
        Vector3 atom1Offset;
        if(l12 < l13) {
            swap(atom1, atom2);
            atom1Offset = atom2Offset;
            r12 = (atom2->position()) - (atom1->position() + atom1Offset);
            r13 = (atom3->position() + atom3Offset) - (atom1->position() + atom1Offset);
        } else {
            swap(atom1, atom3);
            atom1Offset = atom3Offset;
            r12 = (atom2->position() + atom2Offset) - (atom1->position() + atom1Offset);
            r13 = (atom3->position()) - (atom1->position() + atom1Offset);
        }
        dotr12r13 = dot(r12, r13);
        l12Squared = dot(r12, r12);
        l13Squared = dot(r13, r13);
        l12 = sqrt(l12Squared);
        l13 = sqrt(l13Squared);
        angleParam = dotr12r13 / (l12*l13);
        if(angleParam > 1.0 || angleParam < -1.0) {
            // This should never happen, but means that atom2 and atom3 are on top of each other
            return;
        }
        angle = acos(angleParam);
    }

    // Make the potential symmetric, because the neural network might not agree that it is
    if(symmetric && l12 < l13) {
        swap(r12, r13);
        swap(l12Squared, l13Squared);
        swap(l12, l13);
        swap(atom2, atom3);
    }
    double l12MaxOrCutoff = min(l12Max, m_cutoffRadius);
    double l13MaxOrCutoff = min(l13Max, m_cutoffRadius);
    if(l12Squared < l12Min*l12Min || l13Squared < l13Min*l13Min || l12Squared > l12MaxOrCutoff*l12MaxOrCutoff || l13Squared > l13MaxOrCutoff*l13MaxOrCutoff) {
        return;
    }
    // This should not happen if the potential is parameterized between pi/3 and pi
    if(angle < angleMin || angle > angleMax) {
        return;
    }

    double l12Inv = 1 / l12;
    double l13Inv = 1 / l13;

    double shield = 1e-12;
    double invSqrtDotOverLenghtSquared = 1. /
            (
                sqrt(
                    1 - (dotr12r13 * dotr12r13 * (l12Inv * l12Inv * l13Inv * l13Inv)) + shield
                    )
                );

    fann_type input[3];

    input[0] = rescale(l12, l12Min, l12Max);
    input[1] = rescale(l13, l13Min, l13Max);
    input[2] = rescale(angle, angleMin, angleMax);
    double potentialEnergy = rescaleEnergy(fann_run(m_ann, input)[0], energyMin, energyMax);
    potentialEnergy -= energyOffset;


    uint outputIndex = 0;
    FannDerivative::backpropagateDerivative(m_ann, outputIndex);
    double dEdr12 = rescaleEnergyDerivative(m_ann->train_errors[0], l12Min, l12Max, energyMin, energyMax);
    double dEdr13 = rescaleEnergyDerivative(m_ann->train_errors[1], l13Min, l13Max, energyMin, energyMax);
    double dEdangle = rescaleEnergyDerivative(m_ann->train_errors[2], angleMin, angleMax, energyMin, energyMax);

    if(symmetric && fabs(l12 - l13) < 1e-10) {
        // Special case for symmetric potentials
        dEdr13 = dEdr12;
    }

//    cout << atom1->id() << " " << atom2->id() << " " << atom3->id() << endl;
//    cout << setprecision(20);
//    cout << "angle: " << angle << endl;
//    cout << "l12: " << l12 << endl;
//    cout << "l13: " << l13 << endl;
//    cout << "l23: " << Vector3::distanceSquared(atom3->position(), atom2->position()) << endl;
//    cout << "potentialEnergy: " << potentialEnergy << endl;
//    cout << "dEdr12: " << dEdr12 << endl;
//    cout << "dEdr13: " << dEdr13 << endl;
//    cout << "dEdangle: " << dEdangle << endl;

    double dampingFactorDerivativeR12 = 0.0;
    double dampingFactorDerivativeR13 = 0.0;
    double dampingFactorDerivativeAngle = 0.0;
    double dampingFactorR12 = 1.0;
    double dampingFactorR13 = 1.0;
    double dampingFactorAngle = 1.0;
    double totalDampingFactor = 1.0;

    if(damping) {
        double upperLimiter = 0.75;
        double l12DampingMin = l12MaxOrCutoff - upperLimiter;
        double l12DampingMax = l12MaxOrCutoff;
        if(l12 > l12DampingMin && l12 < l12DampingMax) {
            double rij = l12;
            double rd = l12DampingMin;
            double rc = l12DampingMax;
            double exponentialFactor = exp((rij-rd)/(rij-rc));
            dampingFactorR12 *= exponentialFactor * ( (rij-rd)/(rc-rd) + 1 );
            dampingFactorDerivativeR12 = exponentialFactor * ( ( (rij-rd)*(rij + 2*rd - 3*rc) ) / ( (rc - rd)*(rc - rij)*(rc - rij) ) );
//                    cout << "Damping l12!" << endl;
            //        cout << "l12: " << l12 << endl;
            //        cout << "l13: " << l13 << endl;
            //        cout << potentialDampingFactor << endl;
            //        cout << potentialDampingFactorDerivativeR12 << endl;
        }

        double l13DampingMin = l13MaxOrCutoff - upperLimiter;
        double l13DampingMax = l13MaxOrCutoff;
        if(l13 > l13DampingMin && l13 < l13DampingMax) {
            double rij = l13;
            double rd = l13DampingMin;
            double rc = l13DampingMax;
            double exponentialFactor = exp((rij-rd)/(rij-rc));
            dampingFactorR13 *= exponentialFactor * ( (rij-rd)/(rc-rd) + 1 );
            dampingFactorDerivativeR13 = exponentialFactor * ( ( (rij-rd)*(rij + 2*rd - 3*rc) ) / ( (rc - rd)*(rc - rij)*(rc - rij) ) );
//                    cout << "Damping l13!" << endl;
            //        cout << potentialDampingFactor << endl;
            //        cout << potentialDampingFactorDerivativeR13 << endl;
        }

        totalDampingFactor = dampingFactorR12 * dampingFactorR13 * dampingFactorAngle;

        // Note that the derivatives are dependent on the damping factors of the other terms as well
        dEdr12 = potentialEnergy * dampingFactorR13 * dampingFactorAngle * dampingFactorDerivativeR12
                + dEdr12 * totalDampingFactor;
        dEdr13 = potentialEnergy * dampingFactorR12 * dampingFactorAngle * dampingFactorDerivativeR13
                + dEdr13 * totalDampingFactor;
        dEdangle = potentialEnergy * dampingFactorR12 * dampingFactorR13 * dampingFactorDerivativeAngle
                + dEdangle * totalDampingFactor;
    }

    double dangleda = 0;
    double dr12da = 0;
    double dr13da = 0;
    // Atom 1
    for(int a = 0; a < 3; a++) {
        dangleda = -(
                    ((dotr12r13 * r13[a]) * (l12Inv * l13Inv * l13Inv * l13Inv))
                    + ((dotr12r13 * r12[a]) * (l12Inv * l12Inv * l12Inv * l13Inv))
                    + ((-r12[a] - r13[a]) * (l12Inv * l13Inv))
                    ) * invSqrtDotOverLenghtSquared;
        dr12da = -r12[a] * l12Inv;
        dr13da = -r13[a] * l13Inv;

        double forceComp = - dEdr12 * dr12da - dEdr13 * dr13da - dEdangle * dangleda;

        //        forceComp *= softenFactor;

        atom1->addForce(a, forceComp);
    }

    // Atom 2
    for(int a = 0; a < 3; a++) {
        dangleda = -(
                    ((dotr12r13 * (-r12[a])) * (l12Inv * l12Inv * l12Inv * l13Inv))
                    + ((r13[a]) * (l12Inv * l13Inv))
                    ) * invSqrtDotOverLenghtSquared;
        dr12da = r12[a] * l12Inv;
        dr13da = 0;

        double forceComp = - dEdr12 * dr12da - dEdr13 * dr13da - dEdangle * dangleda;

        //        forceComp *= softenFactor;

        atom2->addForce(a, forceComp);

    }

    // Atom 3
    for(int a = 0; a < 3; a++) {
        dangleda = -(
                    ((dotr12r13 * (-r13[a])) * (l12Inv * l13Inv * l13Inv * l13Inv))
                    + ((r12[a]) * (l12Inv * l13Inv))
                    ) * invSqrtDotOverLenghtSquared;
        dr12da = 0;
        dr13da = r13[a] * l13Inv;

        double forceComp = - dEdr12 * dr12da - dEdr13 * dr13da - dEdangle * dangleda;

        //        forceComp *= softenFactor;

        atom3->addForce(a, forceComp);
    }

    atom1->addPotential(totalDampingFactor*potentialEnergy / 3.0);
    atom2->addPotential(totalDampingFactor*potentialEnergy / 3.0);
    atom3->addPotential(totalDampingFactor*potentialEnergy / 3.0);
}
