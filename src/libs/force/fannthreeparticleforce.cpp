#include "fannthreeparticleforce.h"

#include <doublefann.h>
#include <atom.h>
#include <iomanip>

#include <utils/logging.h>

double l12Min = 1.0;
double l12Max = 5.0;
double l13Min = 1.0;
double l13Max = 5.0;
double angleMin = M_PI / 10;
double angleMax = M_PI;
double energyMin = 0.0;
double energyMax = 1.0;

double rescale(double value, double valueMin, double valueMax) {
    return (value - valueMin) / (valueMax - valueMin) * 0.8 + 0.1;
}

double rescaleEnergy(double value, double valueMin, double valueMax) {
    return ((value - 0.1) / 0.8) * (valueMax - valueMin) + valueMin;
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

fann_type* FannThreeParticleForce::testForce(fann_type* input) {
    m_fanntmp = 0;

    // Scaling from ångstrøm to atomic units
    double siToAU = 1.8897;
    double eqDistanceScaled = rescale(1.0*siToAU, l12Min, l12Max);
    double eqAngleScaled = rescale(1.9024, angleMin, angleMax);
    double k1 = 50.0;
    double k2 = 1.0;
    m_fanntmp += k1*(input[0] - eqDistanceScaled) * (input[0] - eqDistanceScaled);
    m_fanntmp += k1*(input[1] - eqDistanceScaled) * (input[1] - eqDistanceScaled);
    m_fanntmp += k2*(input[2] - eqAngleScaled) * (input[2] - eqAngleScaled);
    //    cout << input[2] << " vs. " << eqAngleScaled << endl;
    return &m_fanntmp;
}

void FannThreeParticleForce::calculateAndApplyForce(Atom *atom1, Atom *atom2, Atom *atom3, const Vector3 &atom2Offset, const Vector3 &atom3Offset)
{
    if(!m_ann) {
        warnAboutMissingNetwork();
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

    // Scaling from ångstrøm to atomic units
    double siToAU = 1.8897;

    r12 = r12 * siToAU;
    r13 = r13 * siToAU;

    double shield = 1e-12;

    double dotr12r13 = dot(r12, r13);
    double l12Squared = dot(r12, r12);
    double l13Squared = dot(r13, r13);

    //    if(l12Squared > l12Max*l12Max) {
    //        return;
    //    }
    //    if(l12Squared > l13Max*l13Max) {
    //        return;
    //    }

    double l12 = sqrt(l12Squared);
    double l13 = sqrt(l13Squared);

    // TODO: Use cos angle as parameter instead of angle
    //    double angle2 = acos((l12*l12 + l13*l13 - l23*l23) / (2 * l12 * l13));
    double angleParam = dotr12r13 / (l12*l13);
    if(angleParam > 1.0 || angleParam < -1.0) {
        return; // This should never happen, but means two atom2 and atom3 are on top of each other
    }
    double angle = acos(angleParam);

    if(l12 < l12Min || l12 > l12Max || l13 < l13Min || l13 > l13Max || angle < angleMin || angle > angleMax) {
        return;
    }

    double softenFactor = 1.0;
    double limiter = 0.5;
    if(l12 > l12Max - limiter) {
        softenFactor *= (1 - 1 / limiter * (l12 - (l12Max - limiter)));
    }
    if(l13 > l13Max - limiter) {
        softenFactor *= (1 - 1 / limiter * (l13 - (l13Max - limiter)));
    }
    double angleLimiter = 0.5;
    if(angle > angleMax - angleLimiter) {
        softenFactor *= (1 - 1 / limiter * (angle - (angleMax - angleLimiter)));
    }
//    if(angle < angleMin + angleLimiter) {
//        softenFactor *= (1 - 1 / limiter * (x - (angleMin - angleLimiter)));
//    }

    double l12Inv = 1 / l12;
    double l13Inv = 1 / l13;

    double invSqrtDotOverLenghtSquared = 1. /
            (
                sqrt(
                    1 - (dotr12r13 * dotr12r13 * (l12Inv * l12Inv * l13Inv * l13Inv)) + shield
                    )
                );

    fann_type input[3];
    fann_type *output;

    double energyPlus = 0;
    double energyMinus = 0;

    double h = 1e-8;

    input[0] = rescale(l12, l12Min, l12Max);
    input[1] = rescale(l13, l13Min, l13Max);
    input[2] = rescale(angle, angleMin, angleMax);

    input[0] = rescale(l12 + h, l12Min, l12Max);
    output = fann_run(m_ann, input);
    //    output = testForce(input);
    energyPlus = rescaleEnergy(output[0], energyMin, energyMax);
    input[0] = rescale(l12 - h, l12Min, l12Max);
    output = fann_run(m_ann, input);
    //    output = testForce(input);
    energyMinus = rescaleEnergy(output[0], energyMin, energyMax);
    double dEdr12 = (energyPlus - energyMinus) / (2 * h);

    bool symmetric = true;

    double dEdr13 = 0;
    if(symmetric) {
        input[0] = rescale(l13, l13Min, l13Max);
        input[1] = rescale(l12, l12Min, l12Max);
        input[2] = rescale(angle, angleMin, angleMax);

        input[0] = rescale(l13 + h, l13Min, l13Max);
        output = fann_run(m_ann, input);
        //    output = testForce(input);
        energyPlus = rescaleEnergy(output[0], energyMin, energyMax);
        input[0] = rescale(l13 - h, l13Min, l13Max);
        output = fann_run(m_ann, input);
        //    output = testForce(input);
        energyMinus = rescaleEnergy(output[0], energyMin, energyMax);
        dEdr13 = (energyPlus - energyMinus) / (2 * h);
    } else {
        input[0] = rescale(l12, l12Min, l12Max);
        input[1] = rescale(l13, l13Min, l13Max);
        input[2] = rescale(angle, angleMin, angleMax);

        input[1] = rescale(l13 + h, l13Min, l13Max);
        output = fann_run(m_ann, input);
        //    output = testForce(input);
        energyPlus = rescaleEnergy(output[0], energyMin, energyMax);
        input[1] = rescale(l13 - h, l13Min, l13Max);
        output = fann_run(m_ann, input);
        //    output = testForce(input);
        energyMinus = rescaleEnergy(output[0], energyMin, energyMax);
        dEdr13 = (energyPlus - energyMinus) / (2 * h);
    }

    input[0] = rescale(l12, l12Min, l12Max);
    input[1] = rescale(l13, l13Min, l13Max);
    input[2] = rescale(angle, angleMin, angleMax);

    input[2] = rescale(angle + h, angleMin, angleMax);
    output = fann_run(m_ann, input);
    //    output = testForce(input);
    energyPlus = rescaleEnergy(output[0], energyMin, energyMax);
    input[2] = rescale(angle - h, angleMin, angleMax);
    output = fann_run(m_ann, input);
    //    output = testForce(input);
    energyMinus = rescaleEnergy(output[0], energyMin, energyMax);
    double dEdangle = (energyPlus - energyMinus) / (2 * h);

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

        forceComp *= softenFactor;
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

        forceComp *= softenFactor;
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
        forceComp *= softenFactor;
        atom3->addForce(a, forceComp);

    }

//    double potentialPerAtom = 0.5 * (energyPlus + energyMinus) / 3.0;
    double potentialPerAtom = 100;
    atom1->addPotential(potentialPerAtom);
    atom2->addPotential(potentialPerAtom);
    atom3->addPotential(potentialPerAtom);

}
