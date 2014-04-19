#include "fanntwoparticleforce.h"

#include <doublefann.h>
#include <atom.h>
#include "math/vector3.h"

#include <utils/logging.h>

FannTwoParticleForce::FannTwoParticleForce() :
    m_ann(0),
    m_hasWarnedAboutMissingNetwork(false)
{
}

void FannTwoParticleForce::loadNetwork(const std::string &fileName, const std::string &boundsFilename)
{
    m_ann = fann_create_from_file(fileName.c_str());

    // Read bounds
    ifstream boundsFile(boundsFilename);
    int nInputs;
    boundsFile >> nInputs;
    cout << nInputs << endl;
    int nOutputs;
    boundsFile >> nOutputs;
    cout << nOutputs << endl;

    boundsFile >> m_r12Min;
    boundsFile >> m_r12Max;
    boundsFile >> m_energyMin;
    boundsFile >> m_energyMax;

    cout << "R12 bounds: " << m_r12Min << "," << m_r12Max << endl;
    cout << "Energy bounds: " << m_energyMin << "," << m_energyMax << endl;
}

void FannTwoParticleForce::calculateAndApplyForce(Atom *atom1, Atom *atom2)
{
    calculateAndApplyForce(atom1, atom2, m_zeroVector);
}

bool printed = false;

void FannTwoParticleForce::calculateAndApplyForce(Atom *atom1, Atom *atom2, const Vector3 &atom2Offset)
{
    if(!m_ann) {
        warnAboutMissingNetwork();
        return;
    }

    Vector3 r12 = atom2->position() + atom2Offset - atom1->position();

    double l12Squared = dot(r12, r12);

    // Scaling from ångstrøm to atomic units
    double siToAU = 1.8897;

    l12Squared *= siToAU*siToAU;

    double l12 = sqrt(l12Squared);

    double h = 1e-8;


    double oldl12 = l12;
    bool softForce = false;
    if(l12 > m_r12Max) {
        l12 = m_r12Max;
        softForce = true;
    }

    bool hardForce = false;
    if(l12 < m_r12Min) {
        l12 = m_r12Min;
        hardForce = true;
    }

    fann_type input[1];
    fann_type *output;

    double energyPlus = 0;
    double energyMinus = 0;

    input[0] = rescaleDistance(l12 + h);
    output = fann_run(m_ann, input);
    energyPlus = rescaleEnergy(output[0]);

    input[0] = rescaleDistance(l12 - h);
    output = fann_run(m_ann, input);
    energyMinus = rescaleEnergy(output[0]);

    double dEdr12 = -2.0*(energyPlus - energyMinus) / (2 * h);

    if(hardForce) {
        double diff = oldl12;
        double diff2 = diff*diff;
        double diff4 = diff2*diff2;
        double normal2 = m_r12Min*m_r12Min;
        double normal4 = normal2*normal2;
        dEdr12 += (1 / diff4 - 1 / normal4);
    }

    if(softForce) {
        dEdr12 = 0;
    }

    double dEdr12Normalized = dEdr12 / oldl12;

//    cout << oldl12 << " -> " << dEdr12Normalized << endl;

//    if(softForce) {
//        dEdr12 *= exp(-(l12 - m_r12Max));
//    }

    atom1->addForce(0, -r12.x() * dEdr12Normalized);
    atom1->addForce(1, -r12.y() * dEdr12Normalized);
    atom1->addForce(2, -r12.z() * dEdr12Normalized);

    if(m_isNewtonsThirdLawEnabled) {
        atom2->addForce(0, r12.x() * dEdr12Normalized);
        atom2->addForce(1, r12.y() * dEdr12Normalized);
        atom2->addForce(2, r12.z() * dEdr12Normalized);
        atom2->addPotential((energyPlus + energyMinus) / 4.0);
    }

    atom1->addPotential((energyPlus + energyMinus) / 4.0);
}

void FannTwoParticleForce::warnAboutMissingNetwork()
{
    if(!m_hasWarnedAboutMissingNetwork) {
        LOG(ERROR) << "FANN network not loaded. Cannot apply force." << endl;
        m_hasWarnedAboutMissingNetwork = true;
    }
}

double FannTwoParticleForce::rescaleDistance(double r12)
{
    return (r12 - m_r12Min) / (m_r12Max - m_r12Min) * 0.8 + 0.1;
}

double FannTwoParticleForce::rescaleEnergy(double energy)
{
    return ((energy - 0.1) / 0.8) * (m_energyMax - m_energyMin) + m_energyMin;
}
