#include "fannthreeparticleforce.h"

#include <doublefann.h>
#include <atom.h>
#include <iomanip>

#include <utils/logging.h>
#include <utils/fannderivative.h>

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

    bool symmetric = false;
    bool damping = true;
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
    double siToAU = 1.0;

    r12 = r12 * siToAU;
    r13 = r13 * siToAU;

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
    //    fann_type *output;

    //    double energyPlus = 0;
    //    double energyMinus = 0;

    //    double h = 1e-8;


    //    double l12Input = l12;
    //    double l13Input = l13;
    //    uint l12Index = 0;
    //    uint l13Index = 1;

    input[0] = rescale(l12, l12Min, l12Max);
    input[1] = rescale(l13, l13Min, l13Max);
    input[2] = rescale(angle, angleMin, angleMax);
    double potentialEnergy = rescaleEnergy(fann_run(m_ann, input)[0], energyMin, energyMax);
    potentialEnergy -= energyOffset;
    //    energyPlus = potentialEnergy;
    //    energyMinus = potentialEnergy;


    uint outputIndex = 0;
    FannDerivative::backpropagateDerivative(m_ann, outputIndex);
    double dEdr12 = rescaleEnergyDerivative(m_ann->train_errors[0], l12Min, l12Max, energyMin, energyMax);
    double dEdr13 = rescaleEnergyDerivative(m_ann->train_errors[1], l13Min, l13Max, energyMin, energyMax);
    double dEdangle = rescaleEnergyDerivative(m_ann->train_errors[2], angleMin, angleMax, energyMin, energyMax);
    if(symmetric) {
        input[0] = rescale(l13, l13Min, l13Max);
        input[1] = rescale(l12, l12Min, l12Max);
        input[2] = rescale(angle, angleMin, angleMax);
        double potentialEnergy2 = rescaleEnergy(fann_run(m_ann, input)[0], energyMin, energyMax);
        //        potentialEnergy = 0.5 * ( potentialEnergy + potentialEnergy2);
        FannDerivative::backpropagateDerivative(m_ann, outputIndex);
        dEdr13 = rescaleEnergyDerivative(m_ann->train_errors[0], l13Min, l13Max, energyMin, energyMax);
    }

//    cout << atom1->id() << " " << atom2->id() << " " << atom3->id() << endl;
//    cout << setprecision(20);
//    cout << "angle: " << angle << endl;
//    cout << "l12: " << l12 << endl;
//    cout << "l13: " << l13 << endl;
//    cout << "potentialEnergy: " << potentialEnergy << endl;
//    cout << "dEdr12: " << dEdr12 << endl;
//    cout << "dEdr13: " << dEdr13 << endl;
//    cout << "dEdangle: " << dEdangle << endl;



    //    double softenFactor = 1.0;
    //    //    softenFactor = 0.3;
    //    //     TODO: Consider if this is necessary with 1/3 force for three-particle potentials with equal particles
    //    softenFactor = 1.0 / 3.0; // Because we have equal particles

    double potentialDampingFactorDerivativeR12 = 0.0;
    double potentialDampingFactorDerivativeR13 = 0.0;
    double potentialDampingFactorDerivativeAngle = 0.0;
    double potentialDampingFactor = 1.0;

    if(damping) {
        double limiter = 0.5;
        double l12DampingMin = l12Max - limiter;
        double l12DampingMax = l12Max;
        if(l12 > l12DampingMin && l12 < l12DampingMax) {
            double rij = l12;
            double rd = l12DampingMin;
            double rc = l12DampingMax;
            double exponentialFactor = exp((rij-rd)/(rij-rc));
            potentialDampingFactor *= exponentialFactor * ( (rij-rd)/(rc-rd) + 1 );
            potentialDampingFactorDerivativeR12 = exponentialFactor * ( ( (rij-rd)*(rij + 2*rd - 3*rc) ) / ( (rc - rd)*(rc - rij)*(rc - rij) ) );
                    cout << "Damping l12!" << endl;
            //        cout << "l12: " << l12 << endl;
            //        cout << "l13: " << l13 << endl;
            //        cout << potentialDampingFactor << endl;
            //        cout << potentialDampingFactorDerivativeR12 << endl;
        }

        double l13DampingMin = l13Max - limiter;
        double l13DampingMax = l13Max;
        if(l13 > l13DampingMin && l13 < l13DampingMax) {
            double rij = l13;
            double rd = l13DampingMin;
            double rc = l13DampingMax;
            double exponentialFactor = exp((rij-rd)/(rij-rc));
            potentialDampingFactor *= exponentialFactor * ( (rij-rd)/(rc-rd) + 1 );
            potentialDampingFactorDerivativeR13 = exponentialFactor * ( ( (rij-rd)*(rij + 2*rd - 3*rc) ) / ( (rc - rd)*(rc - rij)*(rc - rij) ) );
                    cout << "Damping l13!" << endl;
            //        cout << potentialDampingFactor << endl;
            //        cout << potentialDampingFactorDerivativeR13 << endl;
        }

//        double angleLimiter = M_PI/6;
        double angleLimiter = 0.001;
        if(angleMax < M_PI - 0.01) { // Avoid damping if max angle is pi
            // Upper angle
            double angleUpperDampingMin = angleMax - angleLimiter;
            double angleUpperDampingMax = angleMax;
            if(angle > angleUpperDampingMin && angle < angleUpperDampingMax) {
                double rij = angle;
                double rd = angleUpperDampingMin;
                double rc = angleUpperDampingMax;
                double exponentialFactor = exp((rij-rd)/(rij-rc));
                potentialDampingFactor *= exponentialFactor * ( (rij-rd)/(rc-rd) + 1 );
                potentialDampingFactorDerivativeAngle = exponentialFactor * ( ( (rij-rd)*(rij + 2*rd - 3*rc) ) / ( (rc - rd)*(rc - rij)*(rc - rij) ) );
                            cout << "Damping angle upper!" << endl;
            }
        }

        // Lower angle
        double angleLowerDampingMin = angleMin + angleLimiter;
        double angleLowerDampingMax = angleMin;
        if(angle > angleLowerDampingMax && angle < angleLowerDampingMin) {
            double rij = angle;
            double rd = angleLowerDampingMin;
            double rc = angleLowerDampingMax;
            double exponentialFactor = exp((rij-rd)/(rij-rc));
            potentialDampingFactor *= exponentialFactor * ( (rij-rd)/(rc-rd) + 1 );
            potentialDampingFactorDerivativeAngle = exponentialFactor * ( ( (rij-rd)*(rij + 2*rd - 3*rc) ) / ( (rc - rd)*(rc - rij)*(rc - rij) ) );
                    cout << "Damping angle lower!" << endl;
            //        cout << "WOOO: " << angleMin << " " << angle << " factor: " << potentialDampingFactor << endl;
            //        cout << "rij: " << rij << endl;
            //        cout << "rd: " << rd << endl;
            //        cout << "rc: " << rc << endl;
            //        cout << "exponentialFactor: " << exponentialFactor << endl;
        }

        dEdr12 = dEdr12 * potentialDampingFactor + potentialEnergy * potentialDampingFactorDerivativeR12;
        dEdr13 = dEdr13 * potentialDampingFactor + potentialEnergy * potentialDampingFactorDerivativeR13;
        dEdangle = dEdangle * potentialDampingFactor + potentialEnergy * potentialDampingFactorDerivativeAngle;
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

    //    double potentialPerAtom = 0.5 * (energyPlus + energyMinus) / 3.0;
    //    double potentialPerAtom = 100;
    //    atom1->addPotential(softenFactor*potentialDampingFactor*potentialEnergy);
    //    atom2->addPotential(softenFactor*potentialDampingFactor*potentialEnergy);
    //    atom3->addPotential(softenFactor*potentialDampingFactor*potentialEnergy);
    atom1->addPotential(potentialDampingFactor*potentialEnergy / 3.0);
    atom2->addPotential(potentialDampingFactor*potentialEnergy / 3.0);
    atom3->addPotential(potentialDampingFactor*potentialEnergy / 3.0);
    //    atom1->addPotential(potentialEnergy / 3.0);
    //    atom2->addPotential(potentialEnergy / 3.0);
    //    atom3->addPotential(potentialEnergy / 3.0);

}
