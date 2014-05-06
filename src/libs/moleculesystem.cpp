#include <moleculesystem.h>

// Local headers
#include <atom.h>
#include <integrator/integrator.h>
#include <moleculesystemcell.h>
#include <force/lennardjonesforce.h>
#include <force/twoparticleforce.h>
#include <integrator/velocityverletintegrator.h>
#include <modifier/modifier.h>
#include <utils/logging.h>

#ifdef USE_MPI
#include <filemanager.h>
#include <processor.h>
#endif
#include <progressreporter.h>
#include <force/threeparticleforce.h>

// System headers
#include <fstream>
#include <iomanip>
#include <sys/stat.h>
#include <algorithm>

using namespace std;

MoleculeSystem::MoleculeSystem() :
    m_nDimensions(3),
    pow3nDimensions(pow(3, m_nDimensions)),
    m_isSaveEnabled(false),
    m_isOutputEnabled(true),
    m_areCellsSetUp(false),
    m_temperature(1.0),
    m_averageDisplacement(0.0),
    m_averageSquareDisplacement(0.0),
    m_pressure(0.0),
    m_step(0),
    m_time(0),
    m_skipInitialize(false),
    #ifdef USE_MPI
    m_processor(new Processor(this)),
    #endif
    m_isCalculatePressureEnabled(true),
    m_isCalculatePotentialEnabled(true),
    m_saveEveryNSteps(1),
    m_nSimulationSteps(0),
    m_isFinalTimeStep(false),
    m_isCreateSymlinkEnabled(false),
    m_twoParticleForce(NULL),
    m_threeParticleForce(NULL),
    m_nAtomsTotal(0)
{
    for(int iDim = 0; iDim < 3; iDim++) {
        m_isPeriodicDimension[iDim] = false;
    }
    m_progressReporter = new ProgressReporter("dragly", "somerun");
    m_integrator = new VelocityVerletIntegrator(this);
#ifdef USE_MPI
    m_fileManager = new FileManager(this);
#endif
    m_cellShiftVectors = zeros(pow3nDimensions, m_nDimensions);
    m_boundaries = zeros(2,m_nDimensions);
}

bool MoleculeSystem::load(string fileName) {
    m_skipInitialize = true;
#ifdef USE_MPI
    return m_fileManager->load(fileName);
#else
    (void)fileName;
#endif
    return true;
}

void MoleculeSystem::clearAtoms()
{
    m_atoms.clear();
}

void MoleculeSystem::addAtomsToCorrectCells(vector<Atom *> &atoms)
{
    for(Atom* atom : atoms) {
        Vector3 position = atom->position();
        for(int iDim = 0; iDim < m_nDimensions; iDim++) {
            if(m_isPeriodicDimension[iDim]) {
                double sideLength = (m_boundaries(1,iDim) - m_boundaries(0,iDim));
                position(iDim) = fmod(position(iDim) + sideLength * 10, sideLength);
            }
        }
        int i = position(0) / m_cellLengths(0);
        int j = position(1) / m_cellLengths(1);
        int k = position(2) / m_cellLengths(2);

        if(i < 0 || i >= m_globalCellsPerDimension(0)) {
            LOG(FATAL) << "Atom ended up outside of system in x direction";
        }
        if(i < 0 || i >= m_globalCellsPerDimension(0)) {
            LOG(FATAL) << "Atom ended up outside of system in y direction";
        }
        if(i < 0 || i >= m_globalCellsPerDimension(0)) {
            LOG(FATAL) << "Atom ended up outside of system in z direction";
        }


        MoleculeSystemCell* cell = m_globalCells.at(cellIndex(i,j,k));
        //        cout << "Put atom into " << cell->indices();
        cell->addAtom(atom);
    }
}

void MoleculeSystem::setStep(uint step)
{
    m_step = step;
}

void MoleculeSystem::deleteAtoms() {
    for(MoleculeSystemCell* cell : m_globalCells) {
        for(Atom* atom : cell->atoms()) {
            delete atom;
        }
        cell->clearAtoms();
    }
    //    m_atoms.clear();
}

const vector<MoleculeSystemCell*>& MoleculeSystem::localCells() const {
#ifdef USE_MPI
    return m_processor->localCells();
#else
    return globalCells();
#endif
}
ThreeParticleForce *MoleculeSystem::threeParticleForce() const
{
    return m_threeParticleForce;
}

void MoleculeSystem::setThreeParticleForce(ThreeParticleForce *threeParticleForce)
{
    m_threeParticleForce = threeParticleForce;
}
TwoParticleForce *MoleculeSystem::twoParticleForce() const
{
    return m_twoParticleForce;
}

void MoleculeSystem::setTwoParticleForce(TwoParticleForce *twoParticleForce)
{
    m_twoParticleForce = twoParticleForce;
}

int MoleculeSystem::cellIndex(int xIndex, int yIndex, int zIndex)
{
    int indexToReturn = zIndex * m_globalCellsPerDimension(0) * m_globalCellsPerDimension(1) + yIndex *  m_globalCellsPerDimension(0) + xIndex;
    if(indexToReturn < 0) {
        LOG(ERROR) << "Got index lower than zero. This should not happen. Arguments were " << xIndex << ", " << yIndex << ", " << zIndex;
    }
    return indexToReturn;
}

void MoleculeSystem::setPeriodicity(bool periodicInX, bool periodicInY, bool periodicInZ)
{
    m_isPeriodicDimension[0] = periodicInX;
    m_isPeriodicDimension[1] = periodicInY;
    m_isPeriodicDimension[2] = periodicInZ;
}

int MoleculeSystem::nAtomsTotal()
{
    return m_nAtomsTotal;
}

const vector<MoleculeSystemCell *> &MoleculeSystem::localAndGhostCells() const
{
#ifdef MD_USE_MPI
    return m_processor->localAndGhostCells();
#else
    return m_globalCells;
#endif
}

void MoleculeSystem::updateStatistics()
{
    // Count the total number of atoms
    m_nAtomsTotal = 0;
    for(MoleculeSystemCell* cell : localCells()) {
        m_nAtomsTotal += cell->atoms().size();
    }
#ifdef USE_MPI
//    LOG(INFO) << world.rank() << " has " << m_nAtomsTotal;
    int nAtomsTotalTmp = 0;
    mpi::all_reduce(world, m_nAtomsTotal, nAtomsTotalTmp, std::plus<int>());
    m_nAtomsTotal = nAtomsTotalTmp;
#endif

    // Calculate total drift
    Vector3 totalDrift;
    totalDrift.zeros();
    for(MoleculeSystemCell* cell : localCells()) {
        for(Atom* atom : cell->atoms()) {
            totalDrift += atom->velocity();
        }
    }

    Vector3 averageDrift;
    averageDrift.zeros();
    if(m_nAtomsTotal > 0) {
        totalDrift / m_nAtomsTotal;
    }

    // Calculate kinetic and potential energy
    m_kineticEnergyTotal = 0;
    m_potentialEnergyTotal = 0;
    for(MoleculeSystemCell* cell : localCells()) {
        for(Atom* atom : cell->atoms()) {
            Vector3 nonDriftVelocity = atom->velocity() - averageDrift;
            m_kineticEnergyTotal += 0.5 * atom->mass() * (dot(nonDriftVelocity, nonDriftVelocity));
            m_potentialEnergyTotal += atom->potential();
        }
    }

#ifdef USE_MPI
    double kineticEnergyTotalTmp = 0.0;
    double potentialEnergyTotalTmp = 0.0;
    mpi::all_reduce(world, m_kineticEnergyTotal, kineticEnergyTotalTmp, std::plus<double>());
    mpi::all_reduce(world, m_potentialEnergyTotal, potentialEnergyTotalTmp, std::plus<double>());
    m_kineticEnergyTotal = kineticEnergyTotalTmp;
    m_potentialEnergyTotal = potentialEnergyTotalTmp;
#endif
    if(m_nAtomsTotal > 0) {
        m_temperature = m_kineticEnergyTotal / (3./2. * m_nAtomsTotal);
    }

    // Calculate diffusion constant
    //    double totalEnergy = (m_potentialEnergyTotal + m_kineticEnergyTotal);
    //    cout << "Atoms: " << nAtomsLocal << " of " << nAtomsTotal << ". Temperature: " << setprecision(5) << m_temperature  << ". Etot: " << totalEnergy << endl;
    m_averageDisplacement = 0;
    m_averageSquareDisplacement = 0;
    for(MoleculeSystemCell* cell : localCells()) {
        for(Atom* atom : cell->atoms()) {
            m_averageSquareDisplacement += (dot(atom->displacement(), atom->displacement()));
            m_averageDisplacement += sqrt(m_averageSquareDisplacement);
        }
    }

    if(m_nAtomsTotal > 0) {
        m_averageSquareDisplacement /= m_nAtomsTotal;
        m_averageDisplacement /= m_nAtomsTotal;
    }

#ifdef USE_MPI
    double averageSquareDisplacementTmp = 0.0;
    double averageDisplacementTmp = 0.0;
    mpi::all_reduce(world, m_averageSquareDisplacement, averageSquareDisplacementTmp, std::plus<double>());
    mpi::all_reduce(world, m_averageDisplacement, averageDisplacementTmp, std::plus<double>());
    m_averageSquareDisplacement = averageSquareDisplacementTmp;
    m_averageDisplacement = averageDisplacementTmp;
#endif
    // Calculate pressure
    rowvec sideLengths = m_boundaries.row(1) - m_boundaries.row(0);
    double volume = (sideLengths(0) * sideLengths(1) * sideLengths(2));
    double density = m_nAtomsTotal / volume;
    double volumeThreeInverse = 1. / (3. * volume);
    m_pressure = 0;
    for(MoleculeSystemCell* cell : localCells()) {
        for(Atom* atom : cell->atoms()) {
            m_pressure += volumeThreeInverse * atom->localPressure();
        }
    }
#ifdef USE_MPI
    double pressureTmp = 0.0;
    mpi::all_reduce(world, m_pressure, pressureTmp, std::plus<double>());
    m_pressure = pressureTmp;
#endif
    m_pressure += density * m_temperature;
    //    cout << "Pressure " << m_pressure << endl;
}

void MoleculeSystem::addAtoms(const vector<Atom *>& atoms)
{
    m_atoms.insert(m_atoms.end(), atoms.begin(), atoms.end());
    obeyBoundaries();
}

const vector<Atom *> &MoleculeSystem::atoms() const
{
    return m_atoms;
}

const vector<MoleculeSystemCell *> &MoleculeSystem::globalCells() const
{
    return m_globalCells;
}

void MoleculeSystem::updateForces()
{
    obeyBoundaries();
    //    refreshCellContents();
#ifdef USE_MPI
    m_processor->communicateAtoms();
    refreshCellContents();
    m_processor->communicateAtoms();
#endif
    refreshCellContents();
    for(MoleculeSystemCell* cell : localAndGhostCells()) {
        for(Atom* atom : cell->atoms()) {
            atom->clearForcePotentialPressure();
            atom->clearNeighborAtoms();
        }
    }

#ifdef USE_MPI
    m_processor->clearForcesInNeighborCells();
#endif

    if(m_twoParticleForce) {
        if(shouldTimeStepBeSaved()) {
            m_twoParticleForce->setCalculatePotentialEnabled(m_isCalculatePotentialEnabled);
            m_twoParticleForce->setCalculatePressureEnabled(m_isCalculatePressureEnabled);
        } else {
            m_twoParticleForce->setCalculatePotentialEnabled(false);
            m_twoParticleForce->setCalculatePressureEnabled(false);
        }
    }

    // Calculate two-particle force and update atoms' lists of neighbor atoms
    for(MoleculeSystemCell* cell : localAndGhostCells()) {
        cell->updateTwoParticleForceAndNeighborAtoms();
    }

    // Use the neighbor lists for each atom to calculate the three-particle force
    if(m_threeParticleForce) {
        for(MoleculeSystemCell* cell : localAndGhostCells()) {
            for(Atom* atom1 : cell->atoms()) {
                int counter = 0;
                for(uint jAtom = 0; jAtom < atom1->neighborAtoms().size(); jAtom++) {
                    std::pair<Atom*,Vector3> atom2 = atom1->neighborAtoms()[jAtom];
                    for(uint kAtom = jAtom + 1; kAtom < atom1->neighborAtoms().size(); kAtom++) {
                        std::pair<Atom*,Vector3> atom3 = atom1->neighborAtoms()[kAtom];
                        const Vector3& offsetVector2 = (atom2.second);
                        const Vector3& offsetVector3 = (atom3.second);
                        m_threeParticleForce->calculateAndApplyForce(atom1, atom2.first, atom3.first, offsetVector2, offsetVector3);
                        counter++;
                    }
                }
            }
        }
    }

#ifdef USE_MPI
    //    m_processor->communicateForces();
#endif
}

MoleculeSystemCell* MoleculeSystem::cell(int i, int j, int k) {
    int iPeriodic = (i + m_globalCellsPerDimension(0)) % m_globalCellsPerDimension(0);
    int jPeriodic = (j + m_globalCellsPerDimension(1)) % m_globalCellsPerDimension(1);
    int kPeriodic = (k + m_globalCellsPerDimension(2)) % m_globalCellsPerDimension(2);
    return m_globalCells.at(iPeriodic + jPeriodic * m_globalCellsPerDimension(0) + kPeriodic * m_globalCellsPerDimension(0) * m_globalCellsPerDimension(1));
}

void MoleculeSystem::simulate()
{
#ifdef USE_MPI
    mpi::timer timer;
    timer.restart();
#endif
    int iStep = 0;
    // Set up integrator if m_skipInitialize is not set (because file was loaded)
    if(!m_skipInitialize) {
        if(isOutputEnabled()) {
            cout << "Initializing integrator" << endl;
        }
        m_integrator->initialize();
    }
    updateStatistics();
    if(isOutputEnabled()) {
        cout << "Starting simulation " << endl;
    }
    for(iStep = iStep; iStep < m_nSimulationSteps; iStep++) {
        m_isFinalTimeStep = (iStep == (m_nSimulationSteps - 1));

        applyModifiers();
        m_integrator->stepForward();
        updateStatistics();

        if(m_processor->rank() == 0 && isOutputEnabledForThisStep()) {
            cout << "Step: " << setw(10) << setfill('0') << m_step
                 << std::scientific << setprecision(16)
                 << " TE: " << m_potentialEnergyTotal + m_kineticEnergyTotal
                 << " PE: " << m_potentialEnergyTotal
                 << " KE: " << m_kineticEnergyTotal
                 << " T: " << m_temperature
                 << " run time: " << setprecision(5) << timer.elapsed()
                 << endl;
        }
#ifdef USE_MPI
        if(shouldTimeStepBeSaved()) {
            m_fileManager->save(m_step);
        }
        if(m_isSaveEnabled) {
            if(m_processor->rank() == 0) {
                m_progressReporter->reportProgress((double)iStep / (double)m_nSimulationSteps);
            }
        }
#endif

        // Finalize step
        m_time += m_integrator->timeStep();
        m_step++;
    }
#ifdef USE_MPI
    if(m_isCreateSymlinkEnabled) {
        m_fileManager->setLatestSymlink(m_step-1);
        if(m_processor->rank() == 0) {
            m_progressReporter->reportProgress(1);
        }
        m_progressReporter->reportProgress(1);
    }
#endif
}

bool MoleculeSystem::shouldTimeStepBeSaved()
{
    if(m_isSaveEnabled && ((m_step % m_saveEveryNSteps) == 0 || m_step == 0)) {
        return true;
    }
    if(m_isSaveEnabled && m_isFinalTimeStep) {
        return true;
    }
    return false;
}

void MoleculeSystem::obeyBoundaries() {
    // Boundary conditions
    for(MoleculeSystemCell* cell : localCells()) {
        for(Atom* atom : cell->atoms()) {
            Vector3 position = atom->position();
            for(int iDim = 0; iDim < m_nDimensions; iDim++) {
                if(m_isPeriodicDimension[iDim]) {
                    double sideLength = (m_boundaries(1,iDim) - m_boundaries(0,iDim));
                    position(iDim) = fmod(position(iDim) + sideLength * 10, sideLength);
                    atom->setPosition(position);
                }
            }
        }
    }
}

void MoleculeSystem::setFileManager(FileManager *fileManager)
{
#ifdef USE_MPI
    m_fileManager = fileManager;
#else
    (void)fileManager;
#endif
}

void MoleculeSystem::addModifier(Modifier *modifier)
{
    m_modifiers.push_back(modifier);
}

void MoleculeSystem::removeModifier(Modifier *modifier)
{
    m_modifiers.erase(std::remove(m_modifiers.begin(), m_modifiers.end(), modifier), m_modifiers.end());
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
    LOG(INFO) << "Setting boundaries to\n" << boundaries;
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

void MoleculeSystem::setupCells(double requestedCellLength) {
//    if(!m_twoParticleForce || m_twoParticleForce->cutoffRadius() <= 0) {
//        LOG(WARNING) << "Two-particle force not set or cutoff radius of two-particle force is <= 0. Defaulting to 3x3x3 cells.";
//        // TODO: Allow one big cell instead...
//        double minimumSideLength = INFINITY;
//        for(int iDim = 0; iDim < m_nDimensions; iDim++) {
//            minimumSideLength = min(minimumSideLength, m_boundaries(1, iDim) - m_boundaries(0,iDim));
//        }
//        cutoffRadius = minimumSideLength;
//    } else {
//        cutoffRadius = m_twoParticleForce->cutoffRadius();
//    }
    LOG(INFO) << "Setting up cells for " << m_atoms.size() << " atoms with cutoff radius: " << requestedCellLength;
    for(uint i = 0; i < m_globalCells.size(); i++) {
        MoleculeSystemCell* cellToDelete = m_globalCells.at(i);
        delete cellToDelete;
    }
    m_globalCells.clear();

    int nCellsTotal = 1;
    m_globalCellsPerDimension = zeros<irowvec>(m_nDimensions);
    m_cellLengths = zeros<rowvec>(m_nDimensions);
    rowvec systemSize = m_boundaries.row(1) - m_boundaries.row(0);
    LOG(INFO) << "System size: " << systemSize;
    for(int iDim = 0; iDim < m_nDimensions; iDim++) {
        double totalLength = m_boundaries(1,iDim) - m_boundaries(0,iDim);
        if(requestedCellLength > totalLength || requestedCellLength == 0.0) {
            requestedCellLength = totalLength;
        }
        m_globalCellsPerDimension(iDim) = totalLength / requestedCellLength;
        m_cellLengths(iDim) = totalLength / m_globalCellsPerDimension(iDim);
        nCellsTotal *= m_globalCellsPerDimension(iDim);
    }
    if(isOutputEnabled()) {
        LOG(INFO) << "Dividing space into a total of " << nCellsTotal << " cells";
        LOG(INFO) << "With geometry " << m_globalCellsPerDimension;
    }

    irowvec indices = zeros<irowvec>(m_nDimensions);

    for(int i = 0; i < nCellsTotal; i++) {
        MoleculeSystemCell* cell = new MoleculeSystemCell(this);
        cell->setID(i);

        mat cellBoundaries = m_boundaries;

        rowvec shiftVector = m_cellLengths % indices;

        cellBoundaries.row(0) = shiftVector;
        cellBoundaries.row(1) = shiftVector + m_cellLengths;

        cell->setBoundaries(cellBoundaries);

        cell->setIndices(indices);

        m_globalCells.push_back(cell);

        indices(0) += 1;
        for(uint iDim = 1; iDim < indices.size(); iDim++) {
            if(indices(iDim - 1) > m_globalCellsPerDimension(iDim - 1) - 1) {
                indices(iDim - 1) = 0;
                indices(iDim) += 1;
            }
        }
    }

//    if(m_globalCells.size() < 27) {
//        LOG(FATAL) << "The number of cells can never be less than 27! Got " << globalCells().size() << " cells.";
//    }

    // Find the neighbor cells
    int nNeighbors;
    for(MoleculeSystemCell *cell1 : globalCells()) {
        nNeighbors = 0;
        for(int i = -1; i <= 1; i++) {
            for(int j = -1; j <= 1; j++) {
                for(int k = -1; k <= 1; k++) {
                    irowvec direction = {i, j, k};
                    irowvec shiftVec = (cell1->indices() + direction);
                    rowvec offsetVec = zeros<rowvec>(3);
                    // Boundaries
                    bool skipNeighbor = false;
                    for(uint l = 0; l < 3; l++) {
                        if(m_isPeriodicDimension[l]) {
                            if(shiftVec(l) >= m_globalCellsPerDimension(l)) {
                                shiftVec(l) -= m_globalCellsPerDimension(l);
                                offsetVec(l) += (m_boundaries(1,l) - m_boundaries(0,l));
                            } else if(shiftVec(l) < 0) {
                                shiftVec(l) += m_globalCellsPerDimension(l);
                                offsetVec(l) -= (m_boundaries(1,l) - m_boundaries(0,l));
                            }
                        } else {
                            // Not periodic boundary conditions and neighbor beyond edge
                            if(shiftVec(l) >= m_globalCellsPerDimension(l) || shiftVec(l) < 0) {
                                skipNeighbor = true;
                            }
                        }
                    }
                    if(skipNeighbor) {
                        continue;
                    }
                    MoleculeSystemCell* cell2 = m_globalCells.at(cellIndex(shiftVec(0), shiftVec(1), shiftVec(2)));
                    cell1->addNeighbor(cell2, offsetVec, direction);
                    nNeighbors++;
                }
            }
        }
    }

    m_areCellsSetUp = true;

    addAtomsToCorrectCells(m_atoms);
#ifdef USE_MPI
    m_processor->setupProcessors();
#endif
}

void MoleculeSystem::refreshCellContents() {
    if(!m_areCellsSetUp) {
        LOG(FATAL) << "Cells must be set up before refreshing their contents!";
    }
    vector<Atom*> allAtoms;
    for(MoleculeSystemCell* cell : m_globalCells) {
        allAtoms.insert(allAtoms.end(), cell->atoms().begin(), cell->atoms().end());
        cell->clearAtoms();
    }
//    LOG(INFO) << "Refreshing cell contents with " << allAtoms.size() << " atoms";
    addAtomsToCorrectCells(allAtoms);

    int atomCounter = 0;
    for(MoleculeSystemCell* cell : m_globalCells) {
        atomCounter += cell->atoms().size();
    }
//    LOG(INFO) << "Now holds " << atomCounter << " atoms in all cells.";
}

void MoleculeSystem::setIntegrator(Integrator *integrator)
{
    m_integrator = integrator;
}

//void MoleculeSystem::addTwoParticleForce(TwoParticleForce *force)
//{
//    m_twoParticleForces.push_back(force);
//}

//void MoleculeSystem::addThreeParticleForce(ThreeParticleForce *force) {
//    m_threeParticleForces.push_back(force);
//}

//const vector<ThreeParticleForce *> &MoleculeSystem::threeParticleForces() const
//{
//    return m_threeParticleForces;
//}

//const vector<TwoParticleForce*>& MoleculeSystem::twoParticleForces() const
//{
//    return m_twoParticleForces;
//}

bool MoleculeSystem::isSaveEnabled() const
{
    return m_isSaveEnabled;
}

void MoleculeSystem::setSaveEnabled(bool enabled)
{
    m_isSaveEnabled = enabled;
}

bool MoleculeSystem::isOutputEnabledForThisStep() const
{
    return (m_isOutputEnabled && ((m_step % m_saveEveryNSteps) == 0 || m_isFinalTimeStep));
}

void MoleculeSystem::save(string fileName)
{
#ifdef USE_MPI
    m_fileManager->setOutFileName(fileName);
    m_fileManager->save(m_step);
#else
    (void)fileName;
#endif
}

void MoleculeSystem::setCreateSymlink(bool enabled)
{
    m_isCreateSymlinkEnabled = enabled;
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

void MoleculeSystem::setParticleTypes(const vector<AtomType>& particleTypes) {
    m_particleTypes = particleTypes;
    for(const AtomType& atomType : particleTypes) {
        m_particleTypesByID[atomType.number()] = atomType;
    }
}
