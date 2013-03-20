#ifndef MOLECULESYSTEM_H
#define MOLECULESYSTEM_H

// Forward declarations
class Atom;
class Integrator;
class MoleculeSystemCell;
class TwoParticleForce;
class FileManager;
class Modifier;
class Processor;

#include <src/math/vector3.h>

// System includes
#include <iostream>
#include <vector>
#include <armadillo>
#include <libconfig.h++>
#include <H5Cpp.h>
#include <H5File.h>

#include <boost/mpi.hpp>
namespace mpi = boost::mpi;

using namespace std;
using namespace arma;
using namespace libconfig;

class MoleculeSystem
{
public:
    MoleculeSystem();

    void loadConfiguration(Config* config);

    void addAtoms(const vector<Atom *> &atoms);
    const vector<Atom*> &atoms() const;
    const vector<MoleculeSystemCell*> &cells() const;
    void updateForces();
    void simulate();

    void setBoundaries(double min, double max);
    void setBoundaries(double xMin, double xMax, double yMin, double yMax, double zMin, double zMax);
    void setBoundaries(mat boundaries);

//    void setupCells();
    void setupCells(double minCutLength);

//    double potentialConstant();

    const mat& cellShiftVectors();

    void refreshCellContents();

//    void setPotentialConstant(double potentialConstant);
    void setIntegrator(Integrator *integrator);
    void setInteratomicForce(TwoParticleForce* force);
    TwoParticleForce* interatomicForce();
    void updateStatistics();

    // I/O
    bool isOutputEnabled() const;
    void setOutputEnabled(bool enabled);
    void setOutFileName(string fileName);
    bool isSaveEnabled() const;
    void setSaveEnabled(bool enabled);
    void obeyBoundaries();
    void setFileManager(FileManager* fileManager);

    // Modifiers
    void addModifier(Modifier* modifier);
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
protected:
    vector<Atom*> m_atoms;
    Integrator *m_integrator;

    double m_potentialConstant;

    mat m_boundaries;

    int m_nDimensions;

    mat m_cellShiftVectors;

    int pow3nDimensions;

    irowvec m_nCells;
    rowvec m_cellLengths;

    vector<MoleculeSystemCell*> m_cells;

    Config *m_config;

    TwoParticleForce *m_interatomicForce;
    FileManager *m_fileManager;

    bool m_isSaveEnabled;
    bool m_isOutputEnabled;

    bool m_areCellsSetUp;
    H5::H5File* hdf5File;

    double m_temperature;

    vector<Modifier*> m_modifiers;
    double m_averageDisplacement;
    double m_averageSquareDisplacement;

    double m_pressure;

    uint m_step;

    double m_time;
    bool m_skipInitialize;
    Processor* m_processor;

    double m_kineticEnergyTotal;
    double m_potentialEnergyTotal;

    mpi::communicator world;
    bool m_isCalculatePressureEnabled;
    bool m_isCalculatePotentialEnabled;

    int m_saveEveryNSteps;
    int m_nSimulationSteps;
};

inline void MoleculeSystem::setNSimulationSteps(int nSteps) {
    m_nSimulationSteps = nSteps;
}

inline int MoleculeSystem::nSimulationSteps() {
    return m_nSimulationSteps;
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
    return m_nCells;
}

inline Processor* MoleculeSystem::processor() {
    return m_processor;
}

inline double MoleculeSystem::kineticEnergyTotal()
{
    return m_kineticEnergyTotal;
}

inline double MoleculeSystem::potentialEnergyTotal()
{
    return m_potentialEnergyTotal;
}

#endif // MOLECULESYSTEM_H
