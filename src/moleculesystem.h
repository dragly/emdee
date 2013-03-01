#ifndef MOLECULESYSTEM_H
#define MOLECULESYSTEM_H

// Forward declarations
class Atom;
class Integrator;
class MoleculeSystemCell;
class TwoParticleForce;
class FileManager;
class Modifier;

#include <src/vector3d.h>

// System includes
#include <iostream>
#include <vector>
#include <armadillo>
#include <libconfig.h++>
#include <H5Cpp.h>
#include <H5File.h>

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
    void simulate(int nSimulationSteps);

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

    void setTime(double currentTime);

    void setBoltzmannConstant(double boltzmannConstant);

    bool load(string fileName);

    void setStep(uint step);
    void deleteAtoms();
    void setAverageSquareDisplacement(double averageSquareDisplacement);
    void setAverageDisplacement(double averageDisplacement);
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
};

#endif // MOLECULESYSTEM_H
