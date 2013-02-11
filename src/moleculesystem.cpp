#include <src/moleculesystem.h>
#include <src/atom.h>
#include <src/molecule.h>
#include <src/integrator/integrator.h>
#include <src/moleculesystemcell.h>
#include <src/interatomicforce.h>
#include <src/integrator/velocityverletintegrator.h>

// System headers
#include <fstream>
#include <iomanip>
#include <sys/stat.h>
#include <algorithm>
#include <boost/lexical_cast.hpp>
#include <H5Cpp.h>
#include <H5File.h>
//#include <hdf5.h>

using namespace std;
using namespace H5;

typedef struct s1_t {
    char atomType[3];
    double positionX;
    double positionY;
    double positionZ;
    double velocityX;
    double velocityY;
    double velocityZ;
    double forceX;
    double forceY;
    double forceZ;
    int cellID;
} s1_t;

MoleculeSystem::MoleculeSystem() :
    m_boltzmannConstant(1.0),
    m_potentialConstant(1.0),
    m_nDimensions(3),
    m_outFileName("/tmp/data*.bin"),
    outFileFormat(XyzFormat),
    pow3nDimensions(pow(3, m_nDimensions)),
    m_unitLength(1),
    m_isSaveEnabled(true),
    m_isOutputEnabled(true),
    m_areCellsSetUp(false)
{
    m_interatomicForce = new InteratomicForce();
    m_integrator = new VelocityVerletIntegrator(this);
    m_cellShiftVectors = zeros(pow3nDimensions, m_nDimensions);
    m_boundaries = zeros(2,m_nDimensions);
}

void MoleculeSystem::loadConfiguration(Config *config)
{
    setUnitLength(config->lookup("units.length"));
    double potentialConstant = config->lookup("system.potentialConstant");
    potentialConstant /= m_unitLength;
    setPotentialConstant(potentialConstant);
    config->lookupValue("simulation.saveFileName", m_outFileName);
    m_isSaveEnabled = config->lookup("simulation.saveEnabled");
}

void MoleculeSystem::load(string fileName) {
    if(isOutputEnabled()) {
        cout << fileName << endl;
    }
}

bool MoleculeSystem::save(int step) {
    if(!m_isSaveEnabled) {
        return true;
    }
    // Refreshing to make sure we don't save any wrong cell IDs
    refreshCellContents();

    if(m_outFileName.find(".xyz") != string::npos) {
        return saveXyz(step);
    } else if(m_outFileName.find(".bin") != string::npos) {
        return saveBinary(step);
    } else if(m_outFileName.find(".h5") != string::npos) {
        return saveHDF5(step);
    }
    return false;
}

bool MoleculeSystem::saveXyz(int step) {
    stringstream outStepName;
    outStepName << setw(6) << setfill('0') << step;
    string outFileNameLocal = m_outFileName;
    size_t starPos = m_outFileName.find("*");
    ofstream outFile;
    if(starPos != string::npos) {
        outFileNameLocal.replace(starPos, 1, outStepName.str());
        outFile.open(outFileNameLocal);
    } else {
        outFile.open(outFileNameLocal, ios_base::app);
    }
    //    cout << "Saving xyz data to " << outFileNameLocal  << endl;

    outFile << m_molecules.size() << endl;
    outFile << "Some nice comment" << endl;

    char line[1000];
    for(MoleculeSystemCell* cell : m_cells) {
        for(Molecule* molecule : cell->molecules()) {
            for(Atom* atom : molecule->atoms()) {
                rowvec position = atom->position() * m_unitLength;
                rowvec velocity = atom->velocity() * m_unitLength;
                rowvec force = atom->force() * m_unitLength;
                sprintf(line, " %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %d\n",
                        position(0), position(1), position(2),
                        velocity(0), velocity(1), velocity(2),
                        force(0), force(1), force(2),
                        cell->id());
                //                toPrint << atom->type().abbreviation
                //                        << " " << position(0) << " " << position(1) << " " << position(2)
                //                        << " " << velocity(0) << " " << velocity(1) << " " << velocity(2)
                //                        << " " << force(0) << " " << force(1) << " " << force(2)
                //                        << " " << cell->id()
                //                        << "\n";
                outFile << atom->type().abbreviation
                        << line;
            }
        }
    }
    outFile.close();
    return true;
}

bool MoleculeSystem::saveBinary(int step) {
    stringstream outStepName;
    outStepName << setw(6) << setfill('0') << step;
    string outFileNameLocal = m_outFileName;
    size_t starPos = m_outFileName.find("*");
    ofstream outFile;
    if(starPos != string::npos) {
        outFileNameLocal.replace(starPos, 1, outStepName.str());
        outFile.open(outFileNameLocal, ios::out | ios::binary);
    } else {
        outFile.open(outFileNameLocal, ios::app | ios::out | ios::binary);
    }
    //    cout << "Saving binary data to " << outFileNameLocal  << endl;

    //    outFile << m_molecules.size() << endl;
    //    outFile << "Some nice comment" << endl;

    //    char line[1000];
    for(MoleculeSystemCell* cell : m_cells) {
        for(Molecule* molecule : cell->molecules()) {
            for(Atom* atom : molecule->atoms()) {
                rowvec position = atom->position() * m_unitLength;
                rowvec velocity = atom->velocity() * m_unitLength;
                rowvec force = atom->force() * m_unitLength;
                //                sprintf(line, " %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %d\n",
                //                        position(0), position(1), position(2),
                //                        velocity(0), velocity(1), velocity(2),
                //                        force(0), force(1), force(2),
                //                        cell->id());
                int id = cell->id();
                char atomType[3];
                sprintf(atomType, "Ar");
                outFile.write(atomType, sizeof(atomType));
                outFile.write((char*)&position(0), sizeof(double));
                outFile.write((char*)&position(1), sizeof(double));
                outFile.write((char*)&position(2), sizeof(double));
                outFile.write((char*)&velocity(0), sizeof(double));
                outFile.write((char*)&velocity(1), sizeof(double));
                outFile.write((char*)&velocity(2), sizeof(double));
                outFile.write((char*)&force(0), sizeof(double));
                outFile.write((char*)&force(1), sizeof(double));
                outFile.write((char*)&force(2), sizeof(double));
                outFile.write((char*)&id, sizeof(int));
            }
        }
    }
    outFile.close();
    return true;
}

void MoleculeSystem::setOutFileName(string fileName)
{
    m_outFileName = fileName;
}

bool MoleculeSystem::saveHDF5(int step) {
    stringstream outStepName;
    outStepName << setw(6) << setfill('0') << step;
    string outFileNameLocal = m_outFileName;
    size_t starPos = m_outFileName.find("*");
    if(starPos != string::npos) {
        outFileNameLocal.replace(starPos, 1, outStepName.str());
    }
    cout << "Writing to " << outFileNameLocal << endl;
    /*
      * Initialize the data
      */
    int  i = 0;
    s1_t* s1 = new s1_t[m_atoms.size()];
    for(MoleculeSystemCell* cell : m_cells) {
        for(Molecule* molecule : cell->molecules()) {
            for(Atom* atom : molecule->atoms()) {
                sprintf(s1[i].atomType , "Ar");
                s1[i].positionX = atom->position()(0);
                s1[i].positionY = atom->position()(1);
                s1[i].positionZ = atom->position()(2);
                s1[i].velocityX = atom->velocity()(0);
                s1[i].velocityY = atom->velocity()(1);
                s1[i].velocityZ = atom->velocity()(2);
                s1[i].forceX = atom->force()(0);
                s1[i].forceY = atom->force()(1);
                s1[i].forceZ = atom->force()(2);
                s1[i].cellID = cell->id();
                i++;
            }
        }
    }

    /*
      * Turn off the auto-printing when failure occurs so that we can
      * handle the errors appropriately
      */
    //        Exception::dontPrint();

    /*
      * Create the data space.
      */
    hsize_t dim[] = {m_atoms.size()};   /* Dataspace dimensions */
    DataSpace space( 1, dim );

    /*
      * Create the file.
      */
    H5File* file = new H5File( outFileNameLocal, H5F_ACC_TRUNC );

    /*
      * Create the memory datatype.
      */
    H5::StrType string_type(H5::PredType::C_S1, 3);
    CompType mtype1( sizeof(s1_t) );
    mtype1.insertMember( "atomType", HOFFSET(s1_t, atomType), string_type);
    mtype1.insertMember( "positionX", HOFFSET(s1_t, positionX), PredType::NATIVE_DOUBLE);
    mtype1.insertMember( "positionY", HOFFSET(s1_t, positionY), PredType::NATIVE_DOUBLE);
    mtype1.insertMember( "positionZ", HOFFSET(s1_t, positionZ), PredType::NATIVE_DOUBLE);
    mtype1.insertMember( "velocityX", HOFFSET(s1_t, velocityX), PredType::NATIVE_DOUBLE);
    mtype1.insertMember( "velocityY", HOFFSET(s1_t, velocityY), PredType::NATIVE_DOUBLE);
    mtype1.insertMember( "velocityZ", HOFFSET(s1_t, velocityZ), PredType::NATIVE_DOUBLE);
    mtype1.insertMember( "forceX", HOFFSET(s1_t, forceX), PredType::NATIVE_DOUBLE);
    mtype1.insertMember( "forceY", HOFFSET(s1_t, forceY), PredType::NATIVE_DOUBLE);
    mtype1.insertMember( "forceZ", HOFFSET(s1_t, forceZ), PredType::NATIVE_DOUBLE);
    mtype1.insertMember( "cellID", HOFFSET(s1_t, cellID), PredType::NATIVE_INT);

    /*
      * Create the dataset.
      */
    DataSet* dataset;
    dataset = new DataSet(file->createDataSet("ArrayOfStructures", mtype1, space));

    /*
      * Write data to the dataset;
      */
    dataset->write( s1, mtype1 );

    /*
      * Release resources
      */
    delete dataset;
    delete file;
    return true;
}

//void MoleculeSystem::setupCells()
//{
//    setupCells(m_potentialConstant * 3);
//    if(m_cells.size() < 27) {
//        cerr << "The number of cells can never be less than 27!" << endl;
//        throw new std::logic_error("The number of cells can never be less than 27!");
//    }
//}

void MoleculeSystem::addMolecules(const vector<Molecule *>& molecule)
{
    m_molecules.insert(m_molecules.end(), molecule.begin(), molecule.end());
}

const vector<Molecule *> &MoleculeSystem::molecules() const
{
    return m_molecules;
}

const vector<MoleculeSystemCell *> &MoleculeSystem::cells() const
{
    return m_cells;
}

void MoleculeSystem::updateForces()
{
    //        for(Molecule* molecule : m_molecules) {
    //            molecule->clearForces();
    //        }
    //        //    double kB = boltzmannConstant;
    //        double eps = m_potentialConstant;
    //        double sigma = 4.5;
    //        Atom* atom1;
    //        Atom* atom2;
    //        rowvec rVec;
    //        rowvec otherPosition;
    //        rowvec shortestVec;
    //        rowvec cellShiftVector;
    //        rowvec cellShiftVectorInUse;
    //        double shortestVecSquaredLength;
    //        double r;
    //        double sigmar;
    //        double factor;
    //        rowvec force;
    //        cout << "Updating forces..." << endl;
    //        int nCalculations = 0;
    //        for(uint i = 0; i < m_atoms.size(); i++) {
    //            for(uint j = 0; j < m_atoms.size(); j++) {
    //                atom1 = m_atoms.at(i);
    //                atom2 = m_atoms.at(j);
    //                if(atom1 == atom2) {
    //                    continue;
    //                }
    //    //            rVec = atom2->absolutePosition() - atom1->absolutePosition();
    //                // Minimum image convention
    //                shortestVecSquaredLength = INFINITY;
    //                for(uint iShiftVec = 0; iShiftVec < m_cellShiftVectors.n_rows; iShiftVec++) {
    //                    cellShiftVector = m_cellShiftVectors.row(iShiftVec);
    //                    otherPosition = atom2->absolutePosition() + cellShiftVector;
    //                    rVec = otherPosition - atom1->absolutePosition();
    //                    double rVecSquaredLength = dot(rVec, rVec);
    //                    if(rVecSquaredLength < shortestVecSquaredLength) {
    //                        shortestVecSquaredLength = rVecSquaredLength;
    //                        shortestVec = rVec;
    //                        cellShiftVectorInUse = cellShiftVector;
    //                    }
    //                }
    //                // Check distances to the nearby cells
    //                r = norm(shortestVec, 2);
    //                sigmar = sigma/r;
    //                // TODO Verify force term
    //                factor = - ((24 * eps) / (r*r)) * (2 * pow((sigmar), 12) - pow((sigmar), 6));

    //                force = factor * rVec;
    //                atom1->addForce(force);
    ////                atom2->addForce(-force);
    //                nCalculations++;
    //            }
    //        }
    //        cout << "nCalculations" << endl;
    //        cout << nCalculations << endl;

    refreshCellContents();
    for(Molecule* molecule : m_molecules) {
        molecule->clearForces();
    }
    for(MoleculeSystemCell* cell : m_cells) {
        cell->clearAlreadyCalculatedNeighbors();
    }
    for(MoleculeSystemCell* cell : m_cells) {
        cell->updateForces();
    }
}

void MoleculeSystem::simulate(int nSimulationSteps)
{
    if(!m_areCellsSetUp) {
        cerr << "Cells must be set up before simulating!" << endl;
        throw(new exception());
    }
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
    m_integrator->initialize();
    if(isOutputEnabled()) {
        cout << "Starting simulation " << endl;
    }
    for(int iStep = 1; iStep < nSimulationSteps; iStep++) {
        if(isOutputEnabled()) {
            cout << iStep << ".." << endl;
        }
        m_integrator->stepForward();
        //        cout << "Integrator done" << endl;

        // Boundary conditions
        for(Molecule* molecule : m_molecules) {
            rowvec position = molecule->position();
            for(int iDim = 0; iDim < m_nDimensions; iDim++) {
                double sideLength = (m_boundaries(1,iDim) - m_boundaries(0,iDim));
                if(fabs(molecule->position()(iDim)) > (m_boundaries(1,iDim) + sideLength)
                        || fabs(molecule->position()(iDim)) < (m_boundaries(0,iDim) - sideLength)) {
                    cerr << "Wow! A molecule ended up  outside 2 x boundaries! The time step must be too big. No reason to continue..." << endl;
                    throw(new exception());
                } else if(molecule->position()(iDim) > m_boundaries(1,iDim)) {
                    position(iDim) -= (m_boundaries(1,iDim) - m_boundaries(0,iDim));
                    molecule->setPosition(position);
                } else if(molecule->position()(iDim) < m_boundaries(0,iDim)) {
                    position(iDim) += (m_boundaries(1,iDim) - m_boundaries(0,iDim));
                    molecule->setPosition(position);
                }
            }
        }
        //        updateForces();
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
    setBoundaries(boundaries); // / m_unitLength);
}

void MoleculeSystem::setBoundaries(mat boundaries)
{
    if(isOutputEnabled()) {
        cout << "Setting boundaries" << endl;
    }
    m_boundaries = boundaries;
    irowvec counters = zeros<irowvec>(m_nDimensions);
    for(int i = 0; i < pow3nDimensions; i++) {
        irowvec direction = counters - ones<irowvec>(m_nDimensions);
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
    double minCutLengthUnit = minCutLength; // / m_unitLength;
    for(uint i = 0; m_cells.size(); i++) {
        MoleculeSystemCell* cellToDelete = m_cells.at(i);
        delete cellToDelete;
    }
    m_cells.clear();

    int nCellsTotal = 1;
    m_nCells = zeros<irowvec>(m_nDimensions);
    m_cellLengths = zeros<rowvec>(m_nDimensions);
    for(int iDim = 0; iDim < m_nDimensions; iDim++) {
        double totalLength = m_boundaries(1,iDim) - m_boundaries(0,iDim);
        m_nCells(iDim) = totalLength / minCutLengthUnit;
        m_cellLengths(iDim) = totalLength / m_nCells(iDim);
        nCellsTotal *= m_nCells(iDim);
    }
    if(isOutputEnabled()) {
        cout << "Dividing space into a total of " << nCellsTotal << " cells" << endl;
        cout << "With geometry " << m_nCells << endl;
    }

    irowvec indices = zeros<irowvec>(m_nDimensions);

    vector<int> randomIDs;
    for(int i = 0; i < nCellsTotal; i++) {
        randomIDs.push_back(i);
    }
    random_shuffle(randomIDs.begin(), randomIDs.end());

    for(int i = 0; i < nCellsTotal; i++) {
        MoleculeSystemCell* cell = new MoleculeSystemCell(this);
        cell->setID(randomIDs.at(i));

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
        irowvec counters = zeros<irowvec>(m_nDimensions);
        nNeighbors = 0;
        for(int i = 0; i < pow3nDimensions; i++) {
            irowvec direction = counters - ones<irowvec>(m_nDimensions);
            irowvec shiftVec = (cell1->indices() + direction);
            rowvec offsetVec = zeros<rowvec>(m_nDimensions);
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
            for(int j = 0; j < m_nDimensions; j++) {
                int multiplicator = 1;
                for(int k = j + 1; k < m_nDimensions; k++) {
                    multiplicator *= m_nCells(m_nDimensions - k - 1);
                }
                cellIndex += multiplicator * shiftVec(m_nDimensions - j - 1);
            }
            //            cout << cellIndex << endl;
            //            cout << offsetVec << endl;
            MoleculeSystemCell* cell2 = m_cells.at(cellIndex);
            if(cell2 != cell1 || m_cells.size() == 1) {
                cell1->addNeighbor(cell2, offsetVec);
                nNeighbors++;
            }
            counters(0) += 1;
            for(uint iDim = 1; iDim < indices.size(); iDim++) {
                if(counters(iDim - 1) >= m_nDimensions) {
                    counters(iDim - 1) = 0;
                    counters(iDim) += 1;
                }
            }
        }
    }

    if(m_cells.size() < 27) {
        cerr << "The number of cells can never be less than 27!" << endl;
        throw new std::logic_error("The number of cells can never be less than 27!");
    }

    refreshCellContents();
    m_areCellsSetUp = true;
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
            for(int iDim = 0; iDim < m_nDimensions; iDim++) {
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
    //    m_potentialConstant = m_potentialConstant; // * m_unitLength / unitLength;
    // Finally set the new unitlength
    m_unitLength = unitLength;
}

void MoleculeSystem::setPotentialConstant(double potentialConstant)
{
    m_potentialConstant = potentialConstant; // / m_unitLength;
}

void MoleculeSystem::setIntegrator(Integrator *integrator)
{
    m_integrator = integrator;
}

void MoleculeSystem::setInteratomicForce(InteratomicForce *force)
{
    m_interatomicForce = force;
}

InteratomicForce *MoleculeSystem::interatomicForce()
{
    return m_interatomicForce;
}

bool MoleculeSystem::isSaveEnabled() const
{
    return m_isSaveEnabled;
}

void MoleculeSystem::setSaveEnabled(bool enabled)
{
    m_isSaveEnabled = enabled;
}

bool MoleculeSystem::isOutputEnabled() const
{
    return m_isOutputEnabled;
}

void MoleculeSystem::setOutputEnabled(bool enabled)
{
    m_isOutputEnabled = enabled;
}

double MoleculeSystem::potentialConstant()
{
    return m_potentialConstant;
}

const mat &MoleculeSystem::cellShiftVectors()
{
    return m_cellShiftVectors;
}

