#ifndef MOLECULESYSTEM_H
#define MOLECULESYSTEM_H

// Forward declarations
class Atom;
class Molecule;
class Integrator;
class MoleculeSystemCell;
class InteratomicForce;
class FileManager;
class Modifier;

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

    void addMolecules(const vector<Molecule *> &molecule);
    const vector<Molecule*> &molecules() const;
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
    void setInteratomicForce(InteratomicForce* force);
    InteratomicForce* interatomicForce();
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

    Integrator* integrator() {
        return m_integrator;
    }
    double temperature() {
        return m_temperature;
    }

    void setBoltzmannConstant(double boltzmannConstant);

protected:
    vector<Molecule*> m_molecules;
    vector<Atom*> m_atoms;
//    int nSimulationSteps;
    Integrator *m_integrator;

    double m_boltzmannConstant;
    double m_potentialConstant;

    mat m_boundaries;

    int m_nDimensions;

    mat m_cellShiftVectors;

    int pow3nDimensions;

    irowvec m_nCells;
    rowvec m_cellLengths;

    vector<MoleculeSystemCell*> m_cells;

    Config *m_config;

    InteratomicForce *m_interatomicForce;
    FileManager *m_fileManager;

    bool m_isSaveEnabled;
    bool m_isOutputEnabled;

    bool m_areCellsSetUp;
    H5::H5File* hdf5File;

    double m_temperature;

    vector<Modifier*> m_modifiers;
};

#endif // MOLECULESYSTEM_H
