#ifndef MOLECULESYSTEM_H
#define MOLECULESYSTEM_H

// Forward declarations
class Atom;
class Integrator;
class MoleculeSystemCell;
class TwoParticleForce;
class ThreeParticleForce;
class FileManager;
class Modifier;
class Processor;
class ProgressReporter;
class SingleParticleForce;

#include <math/vector3.h>
#include <atomtype.h>

// System includes
#include <iostream>
#include <vector>
#include <unordered_map>
#include <armadillo>
//#include <H5Cpp.h>
//#include <H5File.h>

#ifdef USE_MPI
#include <libconfig.h++>
#include <boost/mpi.hpp>
namespace mpi = boost::mpi;
using namespace libconfig;
#endif

using namespace std;
using namespace arma;

class MoleculeSystem
{
public:
    MoleculeSystem();

    void addAtoms(const vector<Atom *> &atoms);
    const vector<Atom*> &atoms() const;
    const vector<MoleculeSystemCell*> &globalCells() const;
    const vector<AtomType> &particleTypes() const;
    void updateForces();
    void simulate();

    void setBoundaries(double min, double max);
    void setBoundaries(double xMin, double xMax, double yMin, double yMax, double zMin, double zMax);
    void setBoundaries(mat boundaries);

//    void setupCells();
    void setupCells();

//    double potentialConstant();

    const mat& cellShiftVectors();

    void refreshCellContents();

//    void setPotentialConstant(double potentialConstant);
    void setIntegrator(Integrator *integrator);
//    void addTwoParticleForce(TwoParticleForce* force);
//    const vector<TwoParticleForce*>& twoParticleForces() const;
    void updateStatistics();

    // I/O
    bool isOutputEnabled() const;
    void setOutputEnabled(bool enabled);
//    void setOutFileName(string fileName);
    bool isSaveEnabled() const;
    void setSaveEnabled(bool enabled);
    void obeyBoundaries();
    void setFileManager(FileManager* fileManager);

    // Modifiers
    void addModifier(Modifier* modifier);
    void removeModifier(Modifier* modifier);
    void applyModifiers();

    // Getters (fast)
    Integrator* integrator() const {
        return m_integrator;
    }
    double kineticEnergyTotal();
    double potentialEnergyTotal();
    double temperature() const {
        return m_temperature;
    }
    double pressure() const {
        return m_pressure;
    }
    double averageDisplacement() const {
        return m_averageDisplacement;
    }
    double averageSquareDisplacement() const {
        return m_averageSquareDisplacement;
    }
    const mat& boundaries() const {
        return m_boundaries;
    }
    double time() const {
        return m_time;
    }
    inline const irowvec& nCells() const;

    void setTime(double currentTime);

    void setBoltzmannConstant(double boltzmannConstant);

    bool load(string fileName);
    void clearAtoms();
    void addAtomsToCorrectCells(vector<Atom*>& atoms);

    void setStep(uint step);
    void deleteAtoms();
    void setAverageSquareDisplacement(double averageSquareDisplacement);
    void setAverageDisplacement(double averageDisplacement);

    void setupProcessors();
    MoleculeSystemCell *cell(int i, int j, int k);

    Processor* processor();

    void setCalculatePressureEnabled(bool enable);
    bool isCalculatePressureEnabled();
    void setCalculatePotentialEnabled(bool enable);
    bool isCalculatePotentialEnabled();
    void setSaveEveryNSteps(int nSteps);
    int saveEveryNSteps();
    bool shouldTimeStepBeSaved();
    void setNSimulationSteps(int nSteps);
    int nSimulationSteps();
    void setProgressReporter(ProgressReporter* progressReporter);
    void addSingleParticleForce(SingleParticleForce* force);
    const vector<SingleParticleForce*>& singleParticleForces() const;
    void setParticleTypes(const vector<AtomType> &particleTypes);
    inline const unordered_map<int, AtomType>& particleTypesById();
//    void addThreeParticleForce(ThreeParticleForce *force);
//    const vector<ThreeParticleForce*>& threeParticleForces() const;
    bool isOutputEnabledForThisStep() const;
    void save(string fileName);
    void setCreateSymlink(bool enabled);
    const vector<MoleculeSystemCell *> &localCells() const;

    ThreeParticleForce *threeParticleForce() const;
    void setThreeParticleForce(ThreeParticleForce *threeParticleForce);

    TwoParticleForce *twoParticleForce() const;
    void setTwoParticleForce(TwoParticleForce *twoParticleForce);

    int cellIndex(int xIndex, int yIndex, int zIndex);
    void setPeriodicity(bool periodicInX, bool periodicInY, bool periodicInZ);

    int nAtomsTotal();

    const vector<MoleculeSystemCell *> &localAndGhostCells() const;

protected:
    vector<Atom*> m_atoms;
    Integrator *m_integrator;

    mat m_boundaries;

    int m_nDimensions;

    mat m_cellShiftVectors;

    int pow3nDimensions;

    irowvec m_globalCellsPerDimension;
    rowvec m_cellLengths;

    vector<MoleculeSystemCell*> m_globalCells;

//    Config *m_config;
//    vector<TwoParticleForce*> m_twoParticleForces;
//    vector<ThreeParticleForce*> m_threeParticleForces;

    bool m_isSaveEnabled;
    bool m_isOutputEnabled;

    bool m_areCellsSetUp;
//    H5::H5File* hdf5File;

    double m_temperature;

    vector<Modifier*> m_modifiers;
    double m_averageDisplacement;
    double m_averageSquareDisplacement;

    double m_pressure;

    uint m_step;

    double m_time;
    bool m_skipInitialize;

    double m_kineticEnergyTotal;
    double m_potentialEnergyTotal;
#ifdef USE_MPI
    FileManager *m_fileManager;
    Processor* m_processor;
    mpi::communicator world;
#endif
    bool m_isCalculatePressureEnabled;
    bool m_isCalculatePotentialEnabled;

    int m_saveEveryNSteps;
    int m_nSimulationSteps;
    bool m_isFinalTimeStep;
    ProgressReporter *m_progressReporter;
    vector<SingleParticleForce*> m_singleParticleForces;
    vector<AtomType> m_particleTypes;
    unordered_map<int,AtomType> m_particleTypesByID;
    bool m_isCreateSymlinkEnabled;
    bool m_isPeriodicDimension[3];

    TwoParticleForce *m_twoParticleForce;
    ThreeParticleForce *m_threeParticleForce;

    int m_nAtomsTotal;
};

inline void MoleculeSystem::addSingleParticleForce(SingleParticleForce *force) {
    m_singleParticleForces.push_back(force);
}

inline const vector<SingleParticleForce *> &MoleculeSystem::singleParticleForces() const
{
    return m_singleParticleForces;
}

inline void MoleculeSystem::setNSimulationSteps(int nSteps) {
    m_nSimulationSteps = nSteps;
}

inline int MoleculeSystem::nSimulationSteps() {
    return m_nSimulationSteps;
}

inline void MoleculeSystem::setProgressReporter(ProgressReporter *progressReporter)
{
    m_progressReporter = progressReporter;
}

inline void MoleculeSystem::setCalculatePressureEnabled(bool enable) {
    m_isCalculatePressureEnabled = enable;
}

inline bool MoleculeSystem::isCalculatePressureEnabled() {
    return m_isCalculatePressureEnabled;
}

inline void MoleculeSystem::setCalculatePotentialEnabled(bool enable) {
    m_isCalculatePotentialEnabled = enable;
}

inline bool MoleculeSystem::isCalculatePotentialEnabled() {
    return m_isCalculatePotentialEnabled;
}

inline void MoleculeSystem::setSaveEveryNSteps(int nSteps) {
    m_saveEveryNSteps = nSteps;
}

inline int MoleculeSystem::saveEveryNSteps() {
    return m_saveEveryNSteps;
}

inline const irowvec &MoleculeSystem::nCells() const {
    return m_globalCellsPerDimension;
}

#ifdef USE_MPI
inline Processor* MoleculeSystem::processor() {
    return m_processor;
}
#endif

inline double MoleculeSystem::kineticEnergyTotal()
{
    return m_kineticEnergyTotal;
}

inline double MoleculeSystem::potentialEnergyTotal()
{
    return m_potentialEnergyTotal;
}

inline const vector<AtomType>& MoleculeSystem::particleTypes() const {
    return m_particleTypes;
}

inline const unordered_map<int, AtomType> &MoleculeSystem::particleTypesById()
{
    return m_particleTypesByID;
}

#endif // MOLECULESYSTEM_H
