#include "fanntwoparticleforce.h"

#include <doublefann.h>
#include <atom.h>
#include "math/vector3.h"

#include <utils/logging.h>
#include <utils/fannderivative.h>

FannTwoParticleForce::FannTwoParticleForce() :
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

    cout << "R12 bounds: " << network.r12Min << "," << network.r12Max << endl;
    cout << "Energy bounds: " << network.energyMin << "," << network.energyMax << endl;
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
    double siToAU = 1.8897;
    r12 = siToAU * r12;

    double l12Squared = dot(r12, r12);

    double l12 = sqrt(l12Squared);

    double h = 1e-8;


    double oldl12 = l12;
    bool softForce = false;
    if(l12 > network->r12Max - 1.0) {
        l12 = network->r12Max - 1.0;
        softForce = true;
    }

    bool hardForce = false;
    if(l12 < network->r12Min) {
        l12 = network->r12Min;
        hardForce = true;
    }

    fann_type input[1];
    fann_type *output;

    double potentialEnergy = 0;
    double energyMinus = 0;

    input[0] = network->rescaleDistance(l12 + h);
    output = fann_run(network->ann, input);
    potentialEnergy = network->rescaleEnergy(output[0]);
    FannDerivative::backpropagateDerivative(network->ann, 0);
    double derivative = network->rescaleEnergyDerivative(network->ann->train_errors[0]);

//    input[0] = network->rescaleDistance(l12 - h);
//    output = fann_run(network->ann, input);
//    energyMinus = network->rescaleEnergy(output[0]);

    double dEdr12 = -1.0*(derivative);

    if(hardForce) {
        double diff = oldl12;
        double diff2 = diff*diff;
        double diff4 = diff2*diff2;
        double normal2 = network->r12Min*network->r12Min;
        double normal4 = normal2*normal2;
        dEdr12 += (1 / diff4 - 1 / normal4);
    }

    if(softForce) {
        double power = 5.0;
        dEdr12 *= 1.0 / pow(oldl12 - l12 + 1, power + 1);
//        dEdr12 *= 0.0;
    }

    double dEdr12Normalized = dEdr12 / oldl12;

//    cout << oldl12 << " -> " << dEdr12Normalized << endl;

//    if(softForce) {
//        dEdr12 *= exp(-(l12 - network->r12Max));
//    }

    atom1->addForce(0, -r12.x() * dEdr12Normalized);
    atom1->addForce(1, -r12.y() * dEdr12Normalized);
    atom1->addForce(2, -r12.z() * dEdr12Normalized);

    if(m_isNewtonsThirdLawEnabled) {
        atom2->addForce(0, r12.x() * dEdr12Normalized);
        atom2->addForce(1, r12.y() * dEdr12Normalized);
        atom2->addForce(2, r12.z() * dEdr12Normalized);
        atom2->addPotential(0.5 * (potentialEnergy + energyMinus) / 2.0);
    }
    atom1->addPotential(0.5 * (potentialEnergy + energyMinus) / 2.0);
}

void FannTwoParticleForce::warnAboutMissingNetwork()
{
    if(!m_hasWarnedAboutMissingNetwork) {
        LOG(ERROR) << "FANN network not loaded. Cannot apply force." << endl;
        m_hasWarnedAboutMissingNetwork = true;
    }
}
