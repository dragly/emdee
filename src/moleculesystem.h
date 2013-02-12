#ifndef MOLECULESYSTEM_H
#define MOLECULESYSTEM_H

// Forward declarations
class Atom;
class Molecule;
class Integrator;
class MoleculeSystemCell;
class InteratomicForce;

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

    enum FileFormat {
        XyzFormat,
        HDF5Format
    };

    void loadConfiguration(Config* config);

    void load(string fileName);
    bool save(int step);
    bool saveXyz(int step);
    void addMolecules(const vector<Molecule *> &molecule);
    const vector<Molecule*> &molecules() const;
    const vector<MoleculeSystemCell*> &cells() const;
    void updateForces();
    void simulate(int nSimulationSteps);

    void setBoundaries(double min, double max);
    void setBoundaries(double xMin, double xMax, double yMin, double yMax, double zMin, double zMax);
    void setBoundaries(mat boundaries);

    bool saveHDF5(int step);
//    void setupCells();
    void setupCells(double minCutLength);

    double potentialConstant();

    const mat& cellShiftVectors();

    void refreshCellContents();

    void setPotentialConstant(double potentialConstant);
    void setIntegrator(Integrator *integrator);
    void setInteratomicForce(InteratomicForce* force);
    InteratomicForce* interatomicForce();
    void updateStatistics();

    // Units
    void setUnitLength(double unitLength);
    void setUnitTime(double unitTime);
    void setUnitEnergy(double unitEnergy);

    // I/O
    bool isOutputEnabled() const;
    void setOutputEnabled(bool enabled);
    bool saveBinary(int step);
    void setOutFileName(string fileName);
    bool isSaveEnabled() const;
    void setSaveEnabled(bool enabled);
protected:
    vector<Molecule*> m_molecules;
    vector<Atom*> m_atoms;
//    int nSimulationSteps;
    Integrator *m_integrator;

    double m_boltzmannConstant;
    double m_potentialConstant;

    mat m_boundaries;

    int m_nDimensions;

    string m_outFileName;
    FileFormat outFileFormat;

    mat m_cellShiftVectors;

    int pow3nDimensions;

    irowvec m_nCells;
    rowvec m_cellLengths;

    vector<MoleculeSystemCell*> m_cells;

    Config *m_config;

    InteratomicForce *m_interatomicForce;

    bool m_isSaveEnabled;
    bool m_isOutputEnabled;

    bool m_areCellsSetUp;
    H5::H5File* hdf5File;

    // Units
    double m_unitLength;
    double m_unitTime;
    double m_unitEnergy;
    double m_unitMass;
};

#endif // MOLECULESYSTEM_H
