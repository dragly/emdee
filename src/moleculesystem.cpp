#include "moleculesystem.h"
#include "atom.h"
#include "molecule.h"
#include "integrator.h"
#include "moleculesystemcell.h"

#include <fstream>
#include <iomanip>
#include <sys/stat.h>
#include <hdf5.h>

using namespace std;

MoleculeSystem::MoleculeSystem() :
//    nSimulationSteps(100),
    boltzmannConstant(2.5),
    m_potentialConstant(0.1),
    nDimensions(3),
    outFileName("out/myfile*.xyz"),
    outFileFormat(XyzFormat),
    pow3nDimensions(pow(3, nDimensions)),
    m_unitLength(1)
{
    integrator = new Integrator(this);
    m_cellShiftVectors = zeros(pow3nDimensions, nDimensions);
    m_boundaries = zeros(2,nDimensions);
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
    outStepName << setw(6) << setfill('0') << step;
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
                    << " " << atom->position()(0) << " " << atom->position()(1) << " " << atom->position()(2)
                    << " " << atom->velocity()(0) << " " << atom->velocity()(1) << " " << atom->velocity()(2)
                    << " " << atom->force()(0) << " " << atom->force()(1) << " " << atom->force()(2)
                    << endl;
        }
    }
    outFile.close();
    return true;
}

bool MoleculeSystem::saveHDF5(string filename) {
    (void)filename;
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

void MoleculeSystem::addMolecules(const vector<Molecule *>& molecule)
{
    m_molecules.insert(m_molecules.end(), molecule.begin(), molecule.end());
}

const vector<Molecule *> &MoleculeSystem::molecules() const
{
    return m_molecules;
}

void MoleculeSystem::updateForces()
{
    //    for(Molecule* molecule : m_molecules) {
    //        molecule->clearForces();
    //    }
    //    //    double kB = boltzmannConstant;
    //    double eps = m_potentialConstant;
    //    double sigma = 4.5;
    //    Atom* atom1;
    //    Atom* atom2;
    //    rowvec rVec;
    //    rowvec otherPosition;
    //    rowvec shortestVec;
    //    rowvec cellShiftVector;
    //    rowvec cellShiftVectorInUse;
    //    double shortestVecSquaredLength;
    //    double r;
    //    double sigmar;
    //    double factor;
    //    rowvec force;
    //    cout << "Updating forces..." << endl;
    //    int nCalculations = 0;
    //    for(uint i = 0; i < m_atoms.size(); i++) {
    //        for(uint j = 0; j < m_atoms.size(); j++) {
    //            atom1 = m_atoms.at(i);
    //            atom2 = m_atoms.at(j);
    //            if(atom1 == atom2) {
    //                continue;
    //            }
    ////            rVec = atom2->absolutePosition() - atom1->absolutePosition();
    //            // Minimum image convention
    ////            shortestVecSquaredLength = INFINITY;
    //            for(uint iShiftVec = 0; iShiftVec < m_cellShiftVectors.n_rows; iShiftVec++) {
    //                cellShiftVector = m_cellShiftVectors.row(iShiftVec);
    //                otherPosition = atom2->absolutePosition() + cellShiftVector;
    //                rVec = otherPosition - atom1->absolutePosition();
    ////                double rVecSquaredLength = dot(rVec, rVec);
    ////                if(rVecSquaredLength < shortestVecSquaredLength) {
    ////                    shortestVecSquaredLength = rVecSquaredLength;
    ////                    shortestVec = rVec;
    ////                    cellShiftVectorInUse = cellShiftVector;
    ////                }
    //                // Check distances to the nearby cells
    //                r = norm(rVec, 2);
    //                sigmar = sigma/r;
    //                // TODO Verify force term
    //                factor = - ((24 * eps) / (r*r)) * (2 * pow((sigmar), 12) - pow((sigmar), 6));

    //                force = factor * rVec;
    //                atom1->addForce(force);
    ////                atom2->addForce(-force);
    //                nCalculations++;
    //            }
    //        }
    //    }
    //    cout << "nCalculations" << endl;
    //    cout << nCalculations << endl;

    refreshCellContents();
    for(MoleculeSystemCell* cell : m_cells) {
        cell->updateForces();
    }
}

void MoleculeSystem::simulate(int nSimulationSteps)
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

    save(0);
    // Set up integrator
    integrator->initialize();
    for(int iStep = 1; iStep < nSimulationSteps; iStep++) {
        cout << "Simulation step " << iStep << endl;
        integrator->stepForward();
        //        cout << "Integrator done" << endl;
        // Boundary conditions
        for(Molecule* molecule : m_molecules) {
            rowvec position = molecule->position();
            for(int iDim = 0; iDim < nDimensions; iDim++) {
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

void MoleculeSystem::setBoundaries(mat boundaries)
{
    m_boundaries = boundaries / m_unitLength;
    irowvec counters = zeros<irowvec>(nDimensions);
    for(int i = 0; i < pow3nDimensions; i++) {
        irowvec direction = counters - ones<irowvec>(nDimensions);
        rowvec lengths = m_boundaries.row(1) - m_boundaries.row(0);
        rowvec directionVec = conv_to<rowvec>::from(direction);
        m_cellShiftVectors.row(i) = lengths % directionVec;
        counters(0) += 1;
        for(uint iDim = 1; iDim < counters.size(); iDim++) {
            if(counters(iDim - 1) > 2) {
                counters(iDim - 1) = 0;
                counters(iDim) += 1;
            }
        }
    }
    //    cout << "Cell shift vectors" << endl;
    //    cout << m_cellShiftVectors << endl;
}

void MoleculeSystem::setupCells(double minCutLength) {
    double minCutLengthUnit = minCutLength / m_unitLength;
    for(uint i = 0; m_cells.size(); i++) {
        MoleculeSystemCell* cellToDelete = m_cells.at(i);
        delete cellToDelete;
    }
    m_cells.clear();

    int nCellsTotal = 1;
    m_nCells = zeros<irowvec>(nDimensions);
    m_cellLengths = zeros<rowvec>(nDimensions);
    for(int iDim = 0; iDim < nDimensions; iDim++) {
        double totalLength = m_boundaries(1,iDim) - m_boundaries(0,iDim);
        m_nCells(iDim) = totalLength / minCutLengthUnit;
        m_cellLengths(iDim) = totalLength / m_nCells(iDim);
        nCellsTotal *= m_nCells(iDim);
    }
    cout << "Dividing space into a total of " << nCellsTotal << " cells" << endl;
    cout << "With geometry " << m_nCells << endl;

    irowvec indices = zeros<irowvec>(nDimensions);
    for(int i = 0; i < nCellsTotal; i++) {
        MoleculeSystemCell* cell = new MoleculeSystemCell(this);

        mat cellBoundaries = m_boundaries;

        rowvec shiftVector = m_cellLengths % indices;

        cellBoundaries.row(0) = shiftVector;
        cellBoundaries.row(1) = shiftVector + m_cellLengths;

        //        cout << "cellBoundaries" << endl;
        //        cout << cellBoundaries << endl;

        cell->setBoundaries(cellBoundaries);

        cell->setIndices(indices);

        m_cells.push_back(cell);

        indices(0) += 1;
        for(uint iDim = 1; iDim < indices.size(); iDim++) {
            if(indices(iDim - 1) > m_nCells(iDim - 1) - 1) {
                indices(iDim - 1) = 0;
                indices(iDim) += 1;
            }
        }
    }

    // Find the neighbor cells
    int nNeighbors;
    for(MoleculeSystemCell *cell1 : m_cells) {
        irowvec counters = zeros<irowvec>(nDimensions);
        nNeighbors = 0;
        for(int i = 0; i < pow3nDimensions; i++) {
            irowvec direction = counters - ones<irowvec>(nDimensions);
            irowvec shiftVec = (cell1->indices() + direction);
            rowvec offsetVec = zeros<rowvec>(nDimensions);
            // Boundaries
            for(uint j = 0; j < shiftVec.n_cols; j++) {
                if(shiftVec(j) >= m_nCells(j)) {
                    shiftVec(j) -= m_nCells(j);
                    offsetVec(j) += (m_boundaries(1,j) - m_boundaries(0,j));
                } else if(shiftVec(j) < 0) {
                    shiftVec(j) += m_nCells(j);
                    offsetVec(j) -= (m_boundaries(1,j) - m_boundaries(0,j));
                }
            }
            int cellIndex = 0;
            for(int j = 0; j < nDimensions; j++) {
                int multiplicator = 1;
                for(int k = j + 1; k < nDimensions; k++) {
                    multiplicator *= m_nCells(nDimensions - k - 1);
                }
                cellIndex += multiplicator * shiftVec(nDimensions - j - 1);
            }
            //            cout << cellIndex << endl;
            //            cout << offsetVec << endl;
            MoleculeSystemCell* cell2 = m_cells.at(cellIndex);
            if(cell2 != cell1) {
                cell1->addNeighbor(cell2, offsetVec);
                nNeighbors++;
            }
            counters(0) += 1;
            for(uint iDim = 1; iDim < indices.size(); iDim++) {
                if(counters(iDim - 1) >= nDimensions) {
                    counters(iDim - 1) = 0;
                    counters(iDim) += 1;
                }
            }
        }
    }

    refreshCellContents();
}

void MoleculeSystem::refreshCellContents() {
    // Add contents (molecules) to the cells
    for(MoleculeSystemCell* cell : m_cells) {
        cell->clearMolecules();
        for(Molecule* molecule : m_molecules) {
            urowvec moreThan = (molecule->position() >= cell->boundaries().row(0));
            urowvec lessThan = (molecule->position() < cell->boundaries().row(1));
            //            cout << "less is more" << endl;
            //            cout << moreThan << endl;
            //            cout << lessThan << endl;
            bool isInside = true;
            for(int iDim = 0; iDim < nDimensions; iDim++) {
                if(!moreThan(iDim) || !lessThan(iDim)) {
                    isInside = false;
                    break;
                }
            }
            if(isInside) {
                cell->addMolecule(molecule);
            }
        }
    }
    //    cout << "The last added cell has " << nNeighbors << " neighbors" << endl;
}

void MoleculeSystem::setUnitLength(double unitLength)
{
    // rescale everything that has been previously set based on previous unit length
    m_potentialConstant = m_potentialConstant * m_unitLength / unitLength;
    // Finally set the new unitlength
    m_unitLength = unitLength;
}

void MoleculeSystem::setPotentialConstant(double potentialConstant)
{
    m_potentialConstant = potentialConstant / m_unitLength;
}

double MoleculeSystem::potentialConstant()
{
    return m_potentialConstant;
}

const mat &MoleculeSystem::cellShiftVectors()
{
    return m_cellShiftVectors;
}

