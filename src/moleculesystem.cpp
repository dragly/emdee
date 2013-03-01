#include <src/moleculesystem.h>

// Local headers
#include <src/atom.h>
#include <src/integrator/integrator.h>
#include <src/moleculesystemcell.h>
#include <src/force/lennardjonesforce.h>
#include <src/force/twoparticleforce.h>
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
    m_potentialConstant(1.0),
    m_nDimensions(3),
    pow3nDimensions(pow(3, m_nDimensions)),
    m_isSaveEnabled(true),
    m_isOutputEnabled(true),
    m_areCellsSetUp(false),
    m_temperature(1.0),
    m_averageDisplacement(0.0),
    m_averageSquareDisplacement(0.0),
    m_pressure(0.0),
    m_step(0),
    m_time(0),
    m_skipInitialize(false)
{
    m_interatomicForce = new LennardJonesForce();
    m_integrator = new VelocityVerletIntegrator(this);
    m_fileManager = new FileManager(this);
    m_cellShiftVectors = zeros(pow3nDimensions, m_nDimensions);
    m_boundaries = zeros(2,m_nDimensions);
}

void MoleculeSystem::loadConfiguration(Config *config)
{
    m_isSaveEnabled = config->lookup("simulation.saveEnabled");
}

bool MoleculeSystem::load(string fileName) {
    m_skipInitialize = true;
    return m_fileManager->load(fileName);
}

void MoleculeSystem::setStep(uint step)
{
    m_step = step;
}

void MoleculeSystem::deleteAtoms() {
    for(Atom* atom : m_atoms) {
        delete atom;
    }
    for(MoleculeSystemCell* cell : m_cells) {
        cell->clearAtoms();
    }
    m_atoms.clear();
}

void MoleculeSystem::updateStatistics()
{
    // Calculate kinetic energy
    double totalKineticEnergy = 0;
    for(Atom* atom : m_atoms) {
        totalKineticEnergy += 0.5 * atom->mass() * (atom->velocity() * atom->velocity());
    }
    m_temperature = totalKineticEnergy / (3./2. * m_atoms.size());
    cout << "Temperature: " << setprecision(25) << m_temperature << endl;

    // Calculate potential energy
    double totalPotentialEnergy = 0;
    for(Atom* atom : m_atoms) {
        totalPotentialEnergy += atom->potential();
    }

    // Calculate diffusion constant
    m_averageDisplacement = 0;
    m_averageSquareDisplacement = 0;
    for(Atom* atom : m_atoms) {
        m_averageSquareDisplacement += (atom->displacement() * atom->displacement());
        m_averageDisplacement += sqrt(m_averageSquareDisplacement);
    }
    m_averageSquareDisplacement /= m_atoms.size();
    m_averageDisplacement /= m_atoms.size();

    // Calculate pressure
    rowvec sideLengths = m_boundaries.row(1) - m_boundaries.row(0);
    double volume = (sideLengths(0) * sideLengths(1) * sideLengths(2));
    double density = m_atoms.size() / volume;
    m_pressure = density * m_temperature;
    double volumeThreeInverse = 1. / (3. * volume);
    for(Atom* atom : m_atoms) {
        m_pressure += volumeThreeInverse * atom->localPressure();
    }
    //    cout << "Pressure " << m_pressure << endl;
}

void MoleculeSystem::addAtoms(const vector<Atom *>& atoms)
{
    m_atoms.insert(m_atoms.end(), atoms.begin(), atoms.end());
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
    obeyBoundaries();
    refreshCellContents();
    for(Atom* atom : m_atoms) {
        atom->clearForcePotentialPressure();
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

    int iStep = 0;
    // Set up integrator if m_skipInitialize is not set (because file was loaded)
    if(!m_skipInitialize) {
        cout << "Initializing integrator" << endl;
        m_integrator->initialize();
        updateStatistics();
        if(m_isSaveEnabled) {
            m_fileManager->save(m_step);
        }
        m_time += m_integrator->timeStep();
        iStep++;
        m_step++;
    }
    if(isOutputEnabled()) {
        cout << "Starting simulation " << endl;
    }
    for(iStep = iStep; iStep < nSimulationSteps; iStep++) {
        if(isOutputEnabled()) {
            cout << "Step " << m_step << ".." << endl;
        }

        applyModifiers();
        m_integrator->stepForward();
        updateStatistics();

        if(m_isSaveEnabled) {
            m_fileManager->save(m_step);
        }
        m_time += m_integrator->timeStep();
        m_step++;
    }
}

void MoleculeSystem::obeyBoundaries() {
    // Boundary conditions
    for(Atom* atom : m_atoms) {
        Vector3 position = atom->position();
        for(int iDim = 0; iDim < m_nDimensions; iDim++) {
            double sideLength = (m_boundaries(1,iDim) - m_boundaries(0,iDim));
            if(fabs(atom->position()(iDim)) > (m_boundaries(1,iDim) + sideLength)
                    || fabs(atom->position()(iDim)) < (m_boundaries(0,iDim) - sideLength)) {
                cerr << "Wow! An atom ended up  outside 2 x boundaries! The time step must be too big. No reason to continue..." << endl;
                throw(new exception());
            } else if(atom->position()(iDim) > m_boundaries(1,iDim)) {
                position(iDim) -= sideLength;
                atom->setPosition(position);
                atom->addDisplacement(sideLength, iDim);
            } else if(atom->position()(iDim) < m_boundaries(0,iDim)) {
                position(iDim) += sideLength;
                atom->setPosition(position);
                atom->addDisplacement(-sideLength, iDim);
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
        cout << "Setting boundaries to" << endl << boundaries;
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
        cout << totalLength << endl;
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
        rowvec counters = zeros<rowvec>(3);
        nNeighbors = 0;
        for(int i = 0; i < pow3nDimensions; i++) {
            rowvec direction = counters - ones<rowvec>(3);
            rowvec shiftVec = (cell1->indices() + direction);
            rowvec offsetVec = zeros<rowvec>(3);
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
        cell->clearAtoms();
    }
    for(Atom* atom : m_atoms) {
        int i = atom->position()(0) / m_cellLengths(0);
        int j = atom->position()(1) / m_cellLengths(1);
        int k = atom->position()(2) / m_cellLengths(2);

        int cellID = k * m_nCells(1) * m_nCells(2) + j *  m_nCells(2) + i;

        //        cout << cellID << endl;

        MoleculeSystemCell* cell = m_cells.at(cellID);
        cell->addAtom(atom);
    }
}

void MoleculeSystem::setIntegrator(Integrator *integrator)
{
    m_integrator = integrator;
}

void MoleculeSystem::setInteratomicForce(TwoParticleForce *force)
{
    m_interatomicForce = force;
}

TwoParticleForce *MoleculeSystem::interatomicForce()
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

void MoleculeSystem::setTime(double currentTime)
{
    m_time = currentTime;
}

void MoleculeSystem::setAverageDisplacement(double averageDisplacement) {
    m_averageDisplacement = averageDisplacement;
}

void MoleculeSystem::setAverageSquareDisplacement(double averageSquareDisplacement) {
    m_averageSquareDisplacement = averageSquareDisplacement;
}
