#include <src/moleculesystem.h>

// Local headers
#include <src/atom.h>
#include <src/molecule.h>
#include <src/integrator/integrator.h>
#include <src/moleculesystemcell.h>
#include <src/interatomicforce.h>
#include <src/integrator/velocityverletintegrator.h>
#include <src/filemanager.h>
#include <src/modifier/modifier.h>

// System headers
#include <fstream>
#include <iomanip>
#include <sys/stat.h>
#include <algorithm>

using namespace std;

MoleculeSystem::MoleculeSystem() :
    m_boltzmannConstant(1.0),
    m_potentialConstant(1.0),
    m_nDimensions(3),
    pow3nDimensions(pow(3, m_nDimensions)),
    m_isSaveEnabled(true),
    m_isOutputEnabled(true),
    m_areCellsSetUp(false),
    m_temperature(1.0)
{
    m_interatomicForce = new InteratomicForce();
    m_integrator = new VelocityVerletIntegrator(this);
    m_fileManager = new FileManager(this);
    m_cellShiftVectors = zeros(pow3nDimensions, m_nDimensions);
    m_boundaries = zeros(2,m_nDimensions);
}

void MoleculeSystem::loadConfiguration(Config *config)
{
//    setUnitLength(config->lookup("units.length"));
//    double potentialConstant = config->lookup("system.potentialConstant");
//    potentialConstant /= m_unitLength;
//    setPotentialConstant(potentialConstant);
    m_isSaveEnabled = config->lookup("simulation.saveEnabled");
}



void MoleculeSystem::updateStatistics()
{
    // Calculate kinetic energy
    double totalKineticEnergy = 0;
    for(Atom* atom : m_atoms) {
        totalKineticEnergy += 0.5 * atom->mass() * dot(atom->velocity(), atom->velocity());
    }
    m_temperature = totalKineticEnergy / (3./2. * m_boltzmannConstant * m_atoms.size());
    cout << "Temperature: " << setprecision(25) << m_temperature << endl;

    // Calculate potential energy
    double totalPotentialEnergy = 0;
    for(Atom* atom : m_atoms) {
        totalPotentialEnergy += atom->potential();
    }
    cout << "totaltPotentialEnergy " << setprecision(25) << totalPotentialEnergy << endl;
    cout << "averagePotentialEnergy " << setprecision(25) << totalPotentialEnergy / m_atoms.size() << endl;

    // Calculate diffusion constant
    rowvec averageDisplacement = zeros<rowvec>(3);
    double averageSquareDisplacement = 0;
    for(Molecule* molecule : m_molecules) {
        averageDisplacement += molecule->displacement();
        averageSquareDisplacement += dot(molecule->displacement(), molecule->displacement());
    }
    cout << "averageDisplacement " << setprecision(25) << averageDisplacement;
    cout << "averageSquareDisplacement " << setprecision(25) << averageSquareDisplacement << endl;

    // Calculate pressure
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

const vector<Atom *> &MoleculeSystem::atoms() const
{
    return m_atoms;
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

    obeyBoundaries();
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
//    cout << "Factor is " << m_unitLength / m_unitTime << endl;
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

    // Set up integrator
    m_integrator->initialize();
    updateStatistics();
    if(m_isSaveEnabled) {
        m_fileManager->save(0);
    }
    if(isOutputEnabled()) {
        cout << "Starting simulation " << endl;
    }
    for(int iStep = 1; iStep < nSimulationSteps; iStep++) {
        if(isOutputEnabled()) {
            cout << iStep << ".." << endl;
        }
        m_integrator->stepForward();
        //        cout << "Integrator done" << endl;
//        updateForces();
//        refreshCellContents();
        //        updateForces();
        applyModifiers();
        updateStatistics();
        if(m_isSaveEnabled) {
            m_fileManager->save(iStep);
        }
    }
}

void MoleculeSystem::obeyBoundaries() {
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
                position(iDim) -= sideLength;
                molecule->setPosition(position);
                molecule->addDisplacement(sideLength, iDim);
            } else if(molecule->position()(iDim) < m_boundaries(0,iDim)) {
                position(iDim) += sideLength;
                molecule->setPosition(position);
                molecule->addDisplacement(-sideLength, iDim);
            }
        }
    }
}

void MoleculeSystem::setFileManager(FileManager *fileManager)
{
    m_fileManager = fileManager;
}

void MoleculeSystem::addModifier(Modifier *modifier)
{
    m_modifiers.push_back(modifier);
}

void MoleculeSystem::setBoltzmannConstant(double boltzmannConstant)
{
    m_boltzmannConstant = boltzmannConstant;
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
    }
    for(Molecule* molecule : m_molecules) {
        int i = molecule->position()(0) / m_cellLengths(0);
        int j = molecule->position()(1) / m_cellLengths(1);
        int k = molecule->position()(2) / m_cellLengths(2);

        int cellID = k * m_nCells(1) * m_nCells(2) + j *  m_nCells(2) + i;

//        cout << cellID << endl;

        MoleculeSystemCell* cell = m_cells.at(cellID);
        cell->addMolecule(molecule);

//        urowvec moreThan = (molecule->position() >= cell->boundaries().row(0));
//        urowvec lessThan = (molecule->position() < cell->boundaries().row(1));
//        //            cout << "less is more" << endl;
//        //            cout << moreThan << endl;
//        //            cout << lessThan << endl;
//        bool isInside = true;
//        for(int iDim = 0; iDim < m_nDimensions; iDim++) {
//            if(!moreThan(iDim) || !lessThan(iDim)) {
//                isInside = false;
//                break;
//            }
//        }
//        if(isInside) {
//            cell->addMolecule(molecule);
//        } else {
//            cout << "wrong!" << endl;
//            exit(919);
//        }
    }
    //    cout << "The last added cell has " << nNeighbors << " neighbors" << endl;
}

//void MoleculeSystem::setPotentialConstant(double potentialConstant)
//{
//    m_potentialConstant = potentialConstant; // / m_unitLength;
//}

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

//double MoleculeSystem::potentialConstant()
//{
//    return m_potentialConstant;
//}

const mat &MoleculeSystem::cellShiftVectors()
{
    return m_cellShiftVectors;
}


void MoleculeSystem::applyModifiers()
{
    for(Modifier* modifier : m_modifiers) {
        modifier->apply();
    }
}
