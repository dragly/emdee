#include "fanntwoparticleforce.h"

#include <doublefann.h>
#include <atom.h>
#include "math/vector3.h"

#include <utils/logging.h>
#include <utils/fannderivative.h>

/*!
 * \class FannTwoParticleForce
 * \brief The FannTwoParticleForce class loads a pre-trained FANN network and uses this
 * to calculate two-particle forces between a given set of atoms.
 */

double inherentEnergy = -0.5;
double tailCorrectionTreshold = 3.0;
double tailCorrectionTresholdUpwards = 0.0;
double beta = 3.0;

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

    cout << "Cutoff radius: " << cutoffRadius() << endl;
    cout << "R12 bounds: " << network.r12Min << "," << network.r12Max << endl;

    network.tailCorrectionMin = network.r12Max - tailCorrectionTreshold;
    network.tailCorrectionMax = network.r12Max + tailCorrectionTresholdUpwards;

    network.tailCorrectionMin = min(network.tailCorrectionMin, cutoffRadius() - tailCorrectionTreshold);
    network.tailCorrectionMax = min(network.tailCorrectionMax, cutoffRadius() + tailCorrectionTresholdUpwards);

    network.headCorrectionMax = network.r12Min;
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
    TwoParticleForce::calculateAndApplyForce(atom1, atom2);
}

bool printed = false;

void FannTwoParticleForce::calculateAndApplyForce(Atom *atom1, Atom *atom2, const Vector3 &atom2Offset)
{
    if(m_networks.empty()) {
        return;
    }

    Vector3 r12 = atom2->position() + atom2Offset - atom1->position();
    double l12Squared = dot(r12, r12);
    if(l12Squared > cutoffRadius()*cutoffRadius()) {
        return;
    }

    double l12 = sqrt(l12Squared);

    const FannTwoParticleNetwork* network = 0;
    for(const FannTwoParticleNetwork& networkA : m_networks) {
        if(l12 <= networkA.r12Min || l12 >= networkA.r12Max) {
            continue;
        }
        if((atom1->type() == networkA.atomType1 && atom2->type() == networkA.atomType2)
                || (atom2->type() == networkA.atomType1 && atom1->type() == networkA.atomType2)) {
            network = &networkA;
            break;
        }
    }
    if(!network) {
        return;
    }

    double dEdr12 = 0.0;

    double potentialEnergy = 0;
    if(l12 > network->tailCorrectionMin) {
        // TODO: Remove duplicate code
        fann_type input[2];
        input[0] = 0.0;
        input[1] = 0.0;
        fann_type *output;

        input[0] = network->rescaleDistance(l12);
        output = fann_run(network->ann, input);
        potentialEnergy = network->rescaleEnergy(output[0]);
        FannDerivative::backpropagateDerivative(network->ann, 0);
        dEdr12 = network->rescaleEnergyDerivative(network->ann->train_errors[0]);

//        dEdr12 = network->tailCorrectionMinEnergy * network->tailCorrectionDampingFactorDerivative(l12) + network->tailCorrectionMinDerivative * network->tailCorrectionDampingFactor(l12);
        dEdr12 = potentialEnergy * network->tailCorrectionDampingFactorDerivative(l12) + dEdr12 * network->tailCorrectionDampingFactor(l12);
//        potentialEnergy = network->tailCorrectionDampingFactor(l12) * network->tailCorrectionMinEnergy;
        potentialEnergy = network->tailCorrectionDampingFactor(l12) * potentialEnergy;
//        cout << potentialEnergy << endl;
//        cout << "Tail!" << endl;q
    } else if(l12 < network->headCorrectionMax) {
        dEdr12 = network->headCorrectionForce(l12);
        potentialEnergy = network->headCorrectionEnergy(l12);
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

    if(isNewtonsThirdLawEnabled()) {
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

double FannTwoParticleNetwork::tailCorrectionDampingFactor(double l12) const
{
    // OLD
    double rij = l12;
    double rd = tailCorrectionMin;
    double rc = tailCorrectionMax;
    double exponentialFactor = exp(beta*(rij-rd)/(rij-rc));
    double dampingFactorR12 = exponentialFactor * ( (rij-rd)/(rc-rd) + 1 );
    return dampingFactorR12;
    // END OLD

    // NEW
//    double rij = l12;
//    double rd = tailCorrectionMin;
//    double powFactor = (rij - rd + 1.0);
//    double powFactor2 = powFactor*powFactor;
//    double powFactor6 = powFactor2*powFactor2*powFactor2;
//    double dampingFactorR12 = 1.0 / powFactor6;
//    return dampingFactorR12;
    // END NEW
}

double FannTwoParticleNetwork::tailCorrectionDampingFactorDerivative(double l12) const
{
    // OLD
    double rij = l12;
    double rd = tailCorrectionMin;
    double rc = tailCorrectionMax;
    double exponentialFactor = exp(beta*(rij-rd)/(rij-rc));
    double dampingFactorDerivativeR12 = beta * exponentialFactor * ( ( (rij-rd)*(rij + 2*rd - 3*rc) ) / ( (rc - rd)*(rc - rij)*(rc - rij) ) );
    return dampingFactorDerivativeR12;
    // END OLD

    // NEW
//    double rd = tailCorrectionMin;
//    double rij = l12;
//    double powFactor = (rij - rd + 1.0);
//    double powFactor2 = powFactor*powFactor;
//    double powFactor6 = powFactor2*powFactor2*powFactor2;
//    double dampingFactorDerivativeR12 = -6.0 / (powFactor6 * powFactor);
//    return dampingFactorDerivativeR12;
    // END NEW
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
    return F0 - (1 / x4 - 1 / x0_4);
}

double FannTwoParticleNetwork::headCorrectionEnergy(double l12) const
{
    double F0 = headCorrectionMaxDerivative;
    double x0 = headCorrectionMax;
    double x = l12;
    double x2 = x*x;
    double x3 = x2*x;
    double x0_2 = x0*x0;
    double x0_3 = x0*x0_2;
    double x0_4 = x0_2*x0_2;
//    return F0 + (1 / x4 - 1 / x0_4);
    // Check this for energy conservation
    double deltaU = -F0*(x - x0) - 1 / (3*x3) + x / (x0_4) - 1 / (3*x0_3) - 1 / x0_3;
    return headCorrectionMaxEnergy + deltaU;
}
