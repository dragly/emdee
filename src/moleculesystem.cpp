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
#include <src/processor.h>
#include <src/progressreporter.h>

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
    m_skipInitialize(false),
    m_processor(new Processor(this)),
    m_isCalculatePressureEnabled(true),
    m_isCalculatePotentialEnabled(true),
    m_saveEveryNSteps(1),
    m_nSimulationSteps(0),
    m_isFinalTimeStep(false)
{
    m_progressReporter = new ProgressReporter("dragly", "somerun");
    m_integrator = new VelocityVerletIntegrator(this);
    m_fileManager = new FileManager(this);
    m_cellShiftVectors = zeros(pow3nDimensions, m_nDimensions);
    m_boundaries = zeros(2,m_nDimensions);
}

bool MoleculeSystem::load(string fileName) {
    m_skipInitialize = true;
    return m_fileManager->load(fileName);
}

void MoleculeSystem::clearAtoms()
{
    m_atoms.clear();
}

void MoleculeSystem::addAtomsToCorrectCells(vector<Atom *> &atoms)
{
    for(Atom* atom : atoms) {
        int i = atom->position()(0) / m_cellLengths(0);
        int j = atom->position()(1) / m_cellLengths(1);
        int k = atom->position()(2) / m_cellLengths(2);

        int cellID = k * m_nCells(0) * m_nCells(1) + j *  m_nCells(0) + i;


        MoleculeSystemCell* cell = m_cells.at(cellID);
        //        cout << "Put atom into " << cell->indices();
        cell->addAtom(atom);
    }
}

void MoleculeSystem::setStep(uint step)
{
    m_step = step;
}

void MoleculeSystem::deleteAtoms() {
    for(MoleculeSystemCell* cell : m_cells) {
        for(Atom* atom : cell->atoms()) {
            delete atom;
        }
        cell->clearAtoms();
    }
    //    m_atoms.clear();
}

void MoleculeSystem::updateStatistics()
{
    // Calculate total drift
    Vector3 totalDrift;
    totalDrift.zeros();
    int nAtomsTotal = 0;
    for(MoleculeSystemCell* cell : m_processor->cells()) {
        nAtomsTotal += cell->atoms().size();
        for(Atom* atom : cell->atoms()) {
            totalDrift += atom->velocity();
        }
    }
    int nAtomsLocal = nAtomsTotal;
    mpi::all_reduce(world, nAtomsTotal, nAtomsTotal, std::plus<int>());

    Vector3 averageDrift = totalDrift / nAtomsTotal;

    // Calculate kinetic and potential energy
    m_kineticEnergyTotal = 0;
    m_potentialEnergyTotal = 0;
    for(MoleculeSystemCell* cell : m_processor->cells()) {
        for(Atom* atom : cell->atoms()) {
            Vector3 nonDriftVelocity = atom->velocity() - averageDrift;
            m_kineticEnergyTotal += 0.5 * atom->mass() * (dot(nonDriftVelocity, nonDriftVelocity));
            m_potentialEnergyTotal += atom->potential();
        }
    }
    mpi::all_reduce(world, m_kineticEnergyTotal, m_kineticEnergyTotal, std::plus<double>());
    mpi::all_reduce(world, m_potentialEnergyTotal, m_potentialEnergyTotal, std::plus<double>());
    m_temperature = m_kineticEnergyTotal / (3./2. * nAtomsTotal);
    cout << "Atoms: " << nAtomsLocal << " of " << nAtomsTotal << ". Temperature: " << setprecision(5) << m_temperature  << endl;

    // Calculate diffusion constant
    if(shouldTimeStepBeSaved()) { // only do these calculations if the time step is saved - these are currently not used elsewhere
        m_averageDisplacement = 0;
        m_averageSquareDisplacement = 0;
        for(MoleculeSystemCell* cell : m_processor->cells()) {
            for(Atom* atom : cell->atoms()) {
                m_averageSquareDisplacement += (dot(atom->displacement(), atom->displacement()));
                m_averageDisplacement += sqrt(m_averageSquareDisplacement);
            }
        }
        m_averageSquareDisplacement /= nAtomsTotal;
        m_averageDisplacement /= nAtomsTotal;
        mpi::all_reduce(world, m_averageSquareDisplacement, m_averageSquareDisplacement, std::plus<double>());
        mpi::all_reduce(world, m_averageDisplacement, m_averageDisplacement, std::plus<double>());

        // Calculate pressure
        rowvec sideLengths = m_boundaries.row(1) - m_boundaries.row(0);
        double volume = (sideLengths(0) * sideLengths(1) * sideLengths(2));
        double density = nAtomsTotal / volume;
        double volumeThreeInverse = 1. / (3. * volume);
        m_pressure = 0;
        for(MoleculeSystemCell* cell : m_processor->cells()) {
            for(Atom* atom : cell->atoms()) {
                m_pressure += volumeThreeInverse * atom->localPressure();
            }
        }
        mpi::all_reduce(world, m_pressure, m_pressure, std::plus<double>());
        m_pressure += density * m_temperature;
        //    cout << "Pressure " << m_pressure << endl;
    }
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
    //    refreshCellContents();
    m_processor->communicateAtoms();
    refreshCellContents();
    for(MoleculeSystemCell* cell : m_processor->cells()) {
        for(Atom* atom : cell->atoms()) {
            atom->clearForcePotentialPressure();
        }
    }

    for(TwoParticleForce* twoParticleForce : m_twoParticleForces) {
        if(shouldTimeStepBeSaved()) {
            twoParticleForce->setCalculatePotentialEnabled(m_isCalculatePotentialEnabled);
            twoParticleForce->setCalculatePressureEnabled(m_isCalculatePressureEnabled);
        } else {
            twoParticleForce->setCalculatePotentialEnabled(false);
            twoParticleForce->setCalculatePressureEnabled(false);
        }
    }
    for(MoleculeSystemCell* cell : m_processor->cells()) {
        cell->updateForces();
    }
}

MoleculeSystemCell* MoleculeSystem::cell(int i, int j, int k) {
    int iPeriodic = (i + m_nCells(0)) % m_nCells(0);
    int jPeriodic = (j + m_nCells(1)) % m_nCells(1);
    int kPeriodic = (k + m_nCells(2)) % m_nCells(2);
    return m_cells.at(iPeriodic + jPeriodic * m_nCells(0) + kPeriodic * m_nCells(0) * m_nCells(1));
}

void MoleculeSystem::simulate()
{
    //    cout << "Factor is " << m_unitLength / m_unitTime << endl;

    // Forget about all cells but our own
    for(MoleculeSystemCell* systemCell : m_cells) {
        bool inMyCells = false;
        for(MoleculeSystemCell* cell : m_processor->cells()) {
            if(cell == systemCell) {
                inMyCells = true;
            }
        }
        if(!inMyCells) {
            for(Atom* atom : systemCell->atoms()) {
                delete atom;
            }
            systemCell->clearAtoms();
        }
    }
    m_atoms.clear();
    for(MoleculeSystemCell* cell : m_processor->cells()) {
        addAtoms(cell->atoms());
    }

    mpi::timer timer;
    timer.restart();
    int iStep = 0;
    // Set up integrator if m_skipInitialize is not set (because file was loaded)
    if(!m_skipInitialize) {
        cout << "Initializing integrator" << endl;
        m_integrator->initialize();
        updateStatistics();
        if(m_isSaveEnabled) {
            m_fileManager->save(m_step);
            if(m_processor->rank() == 0) {
                m_progressReporter->reportProgress(0);
            }
        }

        // Finalize step
        m_time += m_integrator->timeStep();
        iStep++;
        m_step++;
    }
    if(isOutputEnabled()) {
        cout << "Starting simulation " << endl;
    }
    for(iStep = iStep; iStep < m_nSimulationSteps; iStep++) {
        if(isOutputEnabled()) {
            cout << "Step " << m_step << ".." << endl;
        }
        m_isFinalTimeStep = (iStep == (m_nSimulationSteps - 1));

        applyModifiers();
        m_integrator->stepForward();
        updateStatistics();

        if(shouldTimeStepBeSaved()) {
            m_fileManager->save(m_step);
        }
        if(m_isSaveEnabled) {
            if(m_processor->rank() == 0) {
                m_progressReporter->reportProgress((double)iStep / (double)m_nSimulationSteps);
            }
        }

        cout << "Time used so far: " << setprecision(5) << timer.elapsed() << endl;
        // Finalize step
        m_time += m_integrator->timeStep();
        m_step++;
    }
    if(m_isSaveEnabled) {
        m_fileManager->setLatestSymlink(m_step-1);
        if(m_processor->rank() == 0) {
            m_progressReporter->reportProgress(1);
        }
    }
}

bool MoleculeSystem::shouldTimeStepBeSaved()
{
    if(m_isSaveEnabled && (m_step % m_saveEveryNSteps) == 0) {
        return true;
    }
    if(m_isSaveEnabled && m_isFinalTimeStep) {
        return true;
    }
    return false;
}

void MoleculeSystem::obeyBoundaries() {
    // Boundary conditions
    for(MoleculeSystemCell* cell : m_processor->cells()) {
        for(Atom* atom : cell->atoms()) {
            //    for(Atom* atom : m_atoms) {
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
        irowvec counters = zeros<irowvec>(3);
        nNeighbors = 0;
        for(int i = 0; i < pow3nDimensions; i++) {
            irowvec direction = counters - ones<irowvec>(3);
            irowvec shiftVec = (cell1->indices() + direction);
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
                cell1->addNeighbor(cell2, offsetVec, direction);
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

    m_areCellsSetUp = true;

    addAtomsToCorrectCells(m_atoms);

    m_processor->setupProcessors();
}

void MoleculeSystem::refreshCellContents() {
    if(!m_areCellsSetUp) {
        cerr << "Cells must be set up before refreshing their contents!" << endl;
        throw(new exception());
    }
    vector<Atom*> allAtoms;
    for(MoleculeSystemCell* cell : m_cells) {
        allAtoms.insert(allAtoms.end(), cell->atoms().begin(), cell->atoms().end());
        cell->clearAtoms();
    }
//    cout << "Refreshing cell contents with " << allAtoms.size() << " atoms" << endl;
    addAtomsToCorrectCells(allAtoms);
}

void MoleculeSystem::setIntegrator(Integrator *integrator)
{
    m_integrator = integrator;
}

void MoleculeSystem::addTwoParticleForce(TwoParticleForce *force)
{
    m_twoParticleForces.push_back(force);
}

const vector<TwoParticleForce*>& MoleculeSystem::twoParticleForces() const
{
    return m_twoParticleForces;
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
