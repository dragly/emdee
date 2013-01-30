#include "moleculesystem.h"
#include "atom.h"
#include "molecule.h"
#include "integrator.h"

#include <fstream>
#include <iomanip>
#include <sys/stat.h>

using namespace std;

MoleculeSystem::MoleculeSystem() :
    nSimulationSteps(100),
    boltzmannConstant(2.5),
    potentialConstant(0.1),
    nDimensions(3)
{
    integrator = new Integrator(this);
}

void MoleculeSystem::load(string fileName) {
    cout << fileName << endl;
}

bool MoleculeSystem::save(string fileName) {
    string threeLetterExtension = fileName.substr(fileName.length() - 4, 4);
    string fourLetterExtension = fileName.substr(fileName.length() - 5, 5);
    if(threeLetterExtension == ".xyz") {
        return saveXyz(fileName);
    }
    return false;
}

bool MoleculeSystem::saveXyz(string fileName) {
    ofstream outFile;
    outFile.open(fileName);

    outFile << m_molecules.size() << endl;
    outFile << "Some nice comment" << endl;
    for(Molecule* molecule : m_molecules) {
        for(Atom* atom : molecule->atoms()) {
            outFile << atom->type().abbreviation
                    << " " << atom->absolutePosition()(0) << " " << atom->absolutePosition()(1) << " " << atom->absolutePosition()(2)
                    << " " << atom->absoluteVelocity()(0) << " " << atom->absoluteVelocity()(1) << " " << atom->absoluteVelocity()(2)
                    << endl;
        }
    }
    return true;
}

void MoleculeSystem::addMolecules(vector<Molecule *> molecule)
{
    m_molecules.insert(m_molecules.end(), molecule.begin(), molecule.end());
}

const vector<Molecule *> &MoleculeSystem::molecules() const
{
    return m_molecules;
}

void MoleculeSystem::updateForces()
{

    for(Molecule* molecule : m_molecules) {
        molecule->clearForces();
    }
    double kB = boltzmannConstant;
    double eps = potentialConstant;
    double sigma = 4.5;
    Atom* atom1;
    Atom* atom2;
    vec3 rVec;
    double r;
    double sigmar;
    double factor;
    vec3 force;
    for(int iOffset = 0; iOffset < 3; iOffset++) {
        for(int jOffset = 0; jOffset < 3; jOffset++) {
            for(int kOffset = 0; kOffset < 3; kOffset++) {
                vec3 offset;
                vec3 xxOffset;
                xxOffset << iOffset << jOffset << kOffset;
                for(int iDim = 0; iDim < nDimensions; iDim++) {
                    offset(iDim) = (- 1 + xxOffset(iDim)) * (m_boundaries(1,iDim) - m_boundaries(0,iDim)) - m_boundaries(0,iDim);
                }
                for(uint i = 0; i < m_atoms.size(); i++) {
                    for(uint j = i + 1; j < m_atoms.size(); j++) {
                        atom1 = m_atoms.at(i);
                        atom2 = m_atoms.at(j);
                        rVec = atom2->absolutePosition() + offset - atom1->absolutePosition();
                        r = norm(rVec, 2);
                        sigmar = sigma/r;
                        factor = ((24 * eps) / (r*r)) * (2 * pow((sigmar), 12) - pow((sigmar), 6));
                        force = factor * rVec;
                        atom1->addForce(-force);
                        atom2->addForce(force);
                    }
                }
            }
        }
    }
}

void MoleculeSystem::simulate()
{
    // Make output directory
    mkdir("out", S_IRWXU|S_IRGRP|S_IXGRP);

    // Set up atoms list
    m_atoms.clear();
    for(Molecule* molecule : m_molecules) {
        for(Atom* atom : molecule->atoms()) {
            m_atoms.push_back(atom);
        }
    }

    // Set up integrator
    integrator->initialize();
    for(int iStep = 0; iStep < nSimulationSteps; iStep++) {
        cout << "Simulation step " << iStep << endl;
        integrator->stepForward();
        // Boundary conditions
        for(Molecule* molecule : m_molecules) {
            vec3 position = molecule->position();
            for(int iDim = 0; iDim < 3; iDim++) {
                if(molecule->position()(iDim) > m_boundaries(1,iDim)) {
                    position(iDim) -= (m_boundaries(1,iDim) - m_boundaries(0,iDim));
                    molecule->setPosition(position);
                } else if(molecule->position()(iDim) < m_boundaries(0,iDim)) {
                    position(iDim) += (m_boundaries(1,iDim) - m_boundaries(0,iDim));
                    molecule->setPosition(position);
                }
            }
        }
        // Save the data
        stringstream outName;
        outName << "out/mytest" << setw(6) << setfill('0') << iStep << ".xyz";
        save(outName.str());
    }
}

void MoleculeSystem::setBoundaries(double min, double max)
{
    setBoundaries(min, max, min, max, min, max);
}

void MoleculeSystem::setBoundaries(double xMin, double xMax, double yMin, double yMax, double zMin, double zMax)
{
    mat boundaries;
    boundaries << xMin << yMin << zMin << endr
               << xMax << yMax << zMax;
    setBoundaries(boundaries);
}

void MoleculeSystem::setBoundaries(mat::fixed<2,3> boundaries)
{
    m_boundaries = boundaries;
}
