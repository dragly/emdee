#include "fannforce.h"

#include <doublefann.h>
#include <atom.h>

#include <utils/logging.h>

double l12Min = 1.0;
double l12Max = 6.0;
double l13Min = 1.0;
double l13Max = 6.0;
double angleMin = M_PI / 3;
double angleMax = M_PI;

double rescale(double value, double valueMin, double valueMax) {
    return (value - valueMin) / (valueMax - valueMin) * 0.8 + 0.1;
}

FannForce::FannForce() :
    m_ann(0),
    m_hasWarnedAboutMissingNetwork(false)
{
}

void FannForce::warnAboutMissingNetwork() {
    if(!m_hasWarnedAboutMissingNetwork) {
        LOG(ERROR) << "FANN network not loaded. Cannot apply force." << endl;
        m_hasWarnedAboutMissingNetwork = true;
    }
}

void FannForce::loadNetwork(std::string fileName)
{
    m_ann = fann_create_from_file(fileName.c_str());
}

void FannForce::calculateAndApplyForce(Atom *atom1, Atom *atom2, Atom *atom3)
{
    calculateAndApplyForce(atom1, atom2, atom3, m_zeroVector, m_zeroVector);
}

fann_type* FannForce::testForce(fann_type* input) {
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

void FannForce::calculateAndApplyForce(Atom *atom1, Atom *atom2, Atom *atom3, const Vector3 &atom2Offset, const Vector3 &atom3Offset)
{
    (void)atom1;
    (void)atom2;
    (void)atom3;
    (void)atom2Offset;
    (void)atom3Offset;
    if(!m_ann) {
        warnAboutMissingNetwork();
        return;
    }

    if(!(atom1->type().number() == 8 && atom2->type().number() == 1 && atom3->type().number() == 1)) {
        return;
    }

    //    fann_type *output;
    //    fann_type input[3];

    Vector3 r12 = atom2->position() - atom1->position();
    Vector3 r13 = atom3->position() - atom1->position();
    Vector3 r23 = atom3->position() - atom2->position();

    double shield = 1e-12;

    double dotr12r13 = dot(r12, r13);
    double l12Squared = dot(r12, r12);
    double l13Squared = dot(r13, r13);
    double l23Squared = dot(r23, r23);

    // Scaling from ångstrøm to atomic units
    double siToAU = 1.8897;

    l12Squared *= siToAU*siToAU;
    l13Squared *= siToAU*siToAU;
    l23Squared *= siToAU*siToAU;

//    if(l12Squared > l12Max*l12Max) {
//        return;
//    }
//    if(l12Squared > l13Max*l13Max) {
//        return;
//    }

    double l12 = sqrt(l12Squared);
    double l13 = sqrt(l13Squared);
    double l23 = sqrt(l23Squared);

    double l12Inv = 1 / l12;
    double l13Inv = 1 / l12;
    double l23Inv = 1 / l23;

    double invSqrtDotOverLenghtSquared = 1. /
            (
                sqrt(
                    1 - (dotr12r13 * dotr12r13 * (l12Inv * l12Inv * l13Inv * l13Inv)) + shield
                    )
                );


    //    // TODO: Use cos angle as parameter instead of angle
    double angle = acos((l12*l12 + l13*l13 - l23*l23) / (2 * l12 * l13));
//    if(angle < angleMin) {
//        return;
//    }

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
    energyPlus = output[0];
    input[0] = rescale(l12 - h, l12Min, l12Max);
    output = fann_run(m_ann, input);
//    output = testForce(input);
    energyMinus = output[0];
    double dEdr12 = (energyPlus - energyMinus) / (2 * h);

    input[0] = rescale(l12, l12Min, l12Max);
    input[1] = rescale(l13, l13Min, l13Max);
    input[2] = rescale(angle, angleMin, angleMax);

    input[1] = rescale(l13 + h, l13Min, l13Max);
    output = fann_run(m_ann, input);
//    output = testForce(input);
    energyPlus = output[0];
    input[1] = rescale(l13 - h, l13Min, l13Max);
    output = fann_run(m_ann, input);
//    output = testForce(input);
    energyMinus = output[0];
    double dEdr13 = (energyPlus - energyMinus) / (2 * h);

    input[0] = rescale(l12, l12Min, l12Max);
    input[1] = rescale(l13, l13Min, l13Max);
    input[2] = rescale(angle, angleMin, angleMax);

    input[2] = rescale(angle + h, angleMin, angleMax);
    output = fann_run(m_ann, input);
//    output = testForce(input);
    energyPlus = output[0];
    input[2] = rescale(angle - h, angleMin, angleMax);
    output = fann_run(m_ann, input);
//    output = testForce(input);
    energyMinus = output[0];
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
        atom1->addForce(a, forceComp);
    }

    // Atom 2
    for(int a = 0; a < 3; a++) {
        dangleda = -(
                    ((dotr12r13 * (-r12[a])) * (l12Inv * l12Inv * l12Inv * l13Inv))
                    + ((r13[a]) * (l12Inv * l13Inv))
                    ) * invSqrtDotOverLenghtSquared;
        dr12da = r12[a] * l12Inv;

        double forceComp = - dEdr12 * dr12da - dEdr13 * dr13da - dEdangle * dangleda;
        atom2->addForce(a, forceComp);
    }

    // Atom 3
    for(int a = 0; a < 3; a++) {
        dangleda = -(
                    ((dotr12r13 * (-r13[a])) * (l12Inv * l13Inv * l13Inv * l13Inv))
                    + ((r12[a]) * (l12Inv * l13Inv))
                    ) * invSqrtDotOverLenghtSquared;
        dr13da = r13[a] * l13Inv;

        double forceComp = - dEdr12 * dr12da - dEdr13 * dr13da - dEdangle * dangleda;
        atom3->addForce(a, forceComp);
    }

    atom1->addPotential(energyPlus / 3.0);
    atom2->addPotential(energyPlus / 3.0);
    atom3->addPotential(energyPlus / 3.0);
}
