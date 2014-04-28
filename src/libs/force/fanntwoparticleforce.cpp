#include "fanntwoparticleforce.h"

#include <doublefann.h>
#include <atom.h>
#include "math/vector3.h"

#include <utils/logging.h>
#include <utils/fannderivative.h>

FannTwoParticleForce::FannTwoParticleForce() :
    TwoParticleForce(),
    m_hasWarnedAboutMissingNetwork(false)
{
}

void FannTwoParticleForce::addNetwork(const AtomType& atomType1, const AtomType& atomType2,
                                      const std::string &fileName, const std::string &boundsFilename)
{
    FannTwoParticleNetwork network;
    network.ann = fann_create_from_file(fileName.c_str());
    network.atomType1 = atomType1;
    network.atomType2 = atomType2;

    // Read bounds
    ifstream boundsFile(boundsFilename);
    int nInputs;
    boundsFile >> nInputs;
    cout << nInputs << endl;
    int nOutputs;
    boundsFile >> nOutputs;
    cout << nOutputs << endl;

    boundsFile >> network.r12Min;
    boundsFile >> network.r12Max;
    boundsFile >> network.energyMin;
    boundsFile >> network.energyMax;

    network.tailCorrectionMin = network.r12Max - tailCorrectionTreshold;
    network.tailCorrectionMax = cutoffRadius();

    network.headCorrectionMax = network.r12Min;

    cout << "R12 bounds: " << network.r12Min << "," << network.r12Max << endl;
    cout << "Energy bounds: " << network.energyMin << "," << network.energyMax << endl;

    // Find energy at edges
    fann_type input[2];
    input[0] = 0.0;
    input[1] = 0.0;
    double *output;

    // headCorrectionMax
    input[0] = network.rescaleDistance(network.headCorrectionMax);
    output = fann_run(network.ann, input);
    network.headCorrectionMaxEnergy = network.rescaleEnergy(output[0]);
    FannDerivative::backpropagateDerivative(network.ann, 0);
    network.headCorrectionMaxDerivative = network.rescaleEnergyDerivative(network.ann->train_errors[0]);

    // tailCorrectionMin
    input[0] = network.rescaleDistance(network.tailCorrectionMin);
    output = fann_run(network.ann, input);
    network.tailCorrectionMinEnergy = network.rescaleEnergy(output[0]);
    FannDerivative::backpropagateDerivative(network.ann, 0);
    network.tailCorrectionMinDerivative = network.rescaleEnergyDerivative(network.ann->train_errors[0]);

    cout << "Tail correction, min: "
         << network.tailCorrectionMin
         << " max: " << network.tailCorrectionMax
         << " derivative: " << network.tailCorrectionMinDerivative << endl;

//    double F0 = network.tailCorrectionMinForce;
//    double x0 = network.tailCorrectionMin;
//    double x1 = network.tailCorrectionMax;
//    double deltaU1 = -F0*(x1-x0) + F0 / (x1-x0) * (0.5 * x1*x1 - x0*x1 + 0.5 * x0*x0);
//    double U1 = network.tailCorrectionMinEnergy + deltaU1;
    network.energyOffset = 0.0;

    m_networks.push_back(network);
}

void FannTwoParticleForce::calculateAndApplyForce(Atom *atom1, Atom *atom2)
{
    calculateAndApplyForce(atom1, atom2, m_zeroVector);
}

bool printed = false;

void FannTwoParticleForce::calculateAndApplyForce(Atom *atom1, Atom *atom2, const Vector3 &atom2Offset)
{
    if(m_networks.empty()) {
        return;
    }
    const FannTwoParticleNetwork* network = 0;
    for(const FannTwoParticleNetwork& networkA : m_networks) {
        if((atom1->type() == networkA.atomType1 && atom2->type() == networkA.atomType2)
                || (atom2->type() == networkA.atomType1 && atom1->type() == networkA.atomType2)) {
            network = &networkA;
        }
    }
    if(!network) {
        return;
    }
    if(!network->ann) {
        warnAboutMissingNetwork();
        return;
    }

    Vector3 r12 = atom2->position() + atom2Offset - atom1->position();

    // Scaling from ångstrøm to atomic units
    double siToAU = 1.0;
    r12 = siToAU * r12;

    double l12Squared = dot(r12, r12);

    double l12 = sqrt(l12Squared);

    if(l12 > network->tailCorrectionMax) {
        return;
    }

    double dEdr12 = 0.0;

    double potentialEnergy = 0;
    if(l12 > network->tailCorrectionMin) {
        dEdr12 = network->tailCorrectionDerivative(l12);
        potentialEnergy = network->tailCorrectionEnergy(l12);
//        cout << potentialEnergy << endl;
//        cout << "Tail!" << endl;
    } else if(l12 < network->headCorrectionMax) {
        dEdr12 = network->headCorrectionForce(l12);
        potentialEnergy = network->headCorrectionEnergy(l12);
//        cout << "Head!" << endl;
    } else {
        fann_type input[2];
        input[0] = 0.0;
        input[1] = 0.0;
        fann_type *output;

        input[0] = network->rescaleDistance(l12);
        output = fann_run(network->ann, input);
        potentialEnergy = network->rescaleEnergy(output[0]);
        FannDerivative::backpropagateDerivative(network->ann, 0);
        dEdr12 = network->rescaleEnergyDerivative(network->ann->train_errors[0]);
    }
    double force = -1.0*(dEdr12);

    potentialEnergy -= network->energyOffset;

    double dEdr12Normalized = force / l12;

    atom1->addForce(0, -r12.x() * dEdr12Normalized);
    atom1->addForce(1, -r12.y() * dEdr12Normalized);
    atom1->addForce(2, -r12.z() * dEdr12Normalized);

    if(m_isNewtonsThirdLawEnabled) {
        atom2->addForce(0, r12.x() * dEdr12Normalized);
        atom2->addForce(1, r12.y() * dEdr12Normalized);
        atom2->addForce(2, r12.z() * dEdr12Normalized);
        atom2->addPotential(potentialEnergy / 2.0);
    }
    atom1->addPotential(potentialEnergy / 2.0);
}

void FannTwoParticleForce::warnAboutMissingNetwork()
{
    if(!m_hasWarnedAboutMissingNetwork) {
        LOG(ERROR) << "FANN network not loaded. Cannot apply force." << endl;
        m_hasWarnedAboutMissingNetwork = true;
    }
}


double FannTwoParticleNetwork::rescaleDistance(double r12) const
{
    return (r12 - r12Min) / (r12Max - r12Min) * 0.8 + 0.1;
}

double FannTwoParticleNetwork::rescaleEnergy(double energy) const
{
    return ((energy - 0.1) / 0.8) * (energyMax - energyMin) + energyMin;
}

double FannTwoParticleNetwork::rescaleEnergyDerivative(double value) const
{
    return (energyMax - energyMin) / (r12Max - r12Min) * value;
}

double FannTwoParticleNetwork::tailCorrectionEnergy(double l12) const
{
//    double F0 = tailCorrectionMinForce;
//    double x0 = tailCorrectionMin;
//    double x1 = tailCorrectionMax;
//    double x = l12;
////        double deltaU = -F0*(x-x0) + F0 / (x1-x0) * (0.5 * x*x - x0*x + 0.5 * x0*x0);
//    double deltaU = F0*(x0-x)*(x0 - 2*x1 + x) / (2*(x0-x1));
//    return tailCorrectionMinEnergy + deltaU;
    double rij = l12;
    double rd = tailCorrectionMin;
    double rc = tailCorrectionMax;
    double exponentialFactor = exp((rij-rd)/(rij-rc));
    double dampingFactorR12 = exponentialFactor * ( (rij-rd)/(rc-rd) + 1 );
//    dampingFactorDerivativeR12 = exponentialFactor * ( ( (rij-rd)*(rij + 2*rd - 3*rc) ) / ( (rc - rd)*(rc - rij)*(rc - rij) ) );
    return tailCorrectionMinEnergy * dampingFactorR12;
}

double FannTwoParticleNetwork::tailCorrectionDerivative(double l12) const
{
//    double F0 = tailCorrectionMinForce;
//    double x0 = tailCorrectionMin;
//    double x1 = tailCorrectionMax;
//    double x = l12;
//    double F = F0 - F0 * (x - x0) / (x1 - x0);
//    return F;
    double rij = l12;
    double rd = tailCorrectionMin;
    double rc = tailCorrectionMax;
    double exponentialFactor = exp((rij-rd)/(rij-rc));
//    double dampingFactorR12 = exponentialFactor * ( (rij-rd)/(rc-rd) + 1 );
    double dampingFactorDerivativeR12 = exponentialFactor * ( ( (rij-rd)*(rij + 2*rd - 3*rc) ) / ( (rc - rd)*(rc - rij)*(rc - rij) ) );
//    cout << dampingFactorR12 << endl;
    // Because we are going to use a constant energy from here and out, we only get a derivative term from
    // pulling this energy down to zero
    return tailCorrectionMinEnergy * dampingFactorDerivativeR12;
}

double FannTwoParticleNetwork::headCorrectionForce(double l12) const
{
    double F0 = headCorrectionMaxDerivative;
    double x0 = headCorrectionMax;
    double x = l12;
    double x2 = x*x;
    double x4 = x2*x2;
    double x0_2 = x0*x0;
    double x0_4 = x0_2*x0_2;
    return F0 + (1 / x4 - 1 / x0_4);
}

double FannTwoParticleNetwork::headCorrectionEnergy(double l12) const
{
    double F0 = headCorrectionMaxDerivative;
    double x0 = headCorrectionMax;
    double x = l12;
    double x2 = x*x;
    double x3 = x2*x;
    double x4 = x2*x2;
    double x0_2 = x0*x0;
    double x0_3 = x0*x0_2;
    double x0_4 = x0_2*x0_2;
//    return F0 + (1 / x4 - 1 / x0_4);
    double deltaU = -F0*(x - x0) + 1 / (3*x3) + x / (x0_4) - 1 / (3*x0_3) - 1 / x0_3;
    return headCorrectionMaxEnergy + deltaU;
}
