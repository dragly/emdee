#include "moleculesystem.h"
#include "atom.h"
#include "molecule.h"
#include "integrator.h"

#include <fstream>
#include <iomanip>
#include <sys/stat.h>
#include <hdf5.h>

using namespace std;

MoleculeSystem::MoleculeSystem() :
    nSimulationSteps(100),
    boltzmannConstant(2.5),
    potentialConstant(0.1),
    nDimensions(3),
    outFileName("out/myfile.xyz.*"),
    outFileFormat(XyzFormat)
{
    integrator = new Integrator(this);
}

void MoleculeSystem::load(string fileName) {
    cout << fileName << endl;
}

bool MoleculeSystem::save(int step) {
//    string threeLetterExtension = outFileName.substr(outFileName.length() - 4, 4);
    //string fourLetterExtension = outFileName.substr(outFileName.length() - 5, 5);
    if(outFileName.find(".xyz") != string::npos) {
        return saveXyz(step);
    }
    return false;
}

bool MoleculeSystem::saveXyz(int step) {
    stringstream outStepName;
    outStepName << setw(4) << setfill('0') << step;
    string outFileNameLocal = outFileName;
    size_t starPos = outFileName.find("*");
    ofstream outFile;
    if(starPos != string::npos) {
        outFileNameLocal.replace(starPos, 1, outStepName.str());
        outFile.open(outFileNameLocal);
    } else {
        outFile.open(outFileNameLocal, ios_base::app);
    }
    cout << "Saving xyz data to " << outFileNameLocal  << endl;

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
    outFile.close();
    return true;
}

bool MoleculeSystem::saveHDF5(string filename) {
//    hid_t       file, filetype, memtype, strtype, space, dset;
//                                            /* Handles */
//    herr_t      status;
//    hsize_t     dims[1] = {DIM0};
//    Molecule    wdata[DIM0],                /* Write buffer */
//                *rdata;                     /* Read buffer */
//    int         ndims,
//                i;

//    /*
//     * Create a new file using the default properties.
//     */
//    file = H5Fcreate (FILE, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

//    /*
//     * Create variable-length string datatype.
//     */
//    strtype = H5Tcopy (H5T_C_S1);
//    status = H5Tset_size (strtype, H5T_VARIABLE);

//    /*
//     * Create the compound datatype for memory.
//     */
//    memtype = H5Tcreate (H5T_COMPOUND, sizeof (Molecule));
//    status = H5Tinsert (memtype, "Mass", HOFFSET (Molecule, mass), H5T_NATIVE_DOUBLE);

//    space = H5Screate_simple (1, dims, NULL);

//    dset = H5Dcreate (file, DATASET, filetype, space, H5P_DEFAULT, H5P_DEFAULT,
//                H5P_DEFAULT);
//    status = H5Dwrite (dset, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, wdata);
//    /*
//     * Close and release resources.
//     */
//    status = H5Dclose (dset);
//    status = H5Sclose (space);
//    status = H5Tclose (filetype);
//    status = H5Fclose (file);
    return false;
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
        save(iStep);
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
