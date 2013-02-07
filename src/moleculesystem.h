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
    void updateForces();
    void simulate(int nSimulationSteps);

    void setBoundaries(double min, double max);
    void setBoundaries(double xMin, double xMax, double yMin, double yMax, double zMin, double zMax);
    void setBoundaries(mat boundaries);

    bool saveHDF5(string filename);
    void setupCells();
    void setupCells(double minCutLength);

    double potentialConstant();

    const mat& cellShiftVectors();

    void refreshCellContents();

    void setUnitLength(double unitLength);
    void setPotentialConstant(double potentialConstant);
    void setIntegrator(Integrator *integrator);
    void setInteratomicForce(InteratomicForce* force);
    InteratomicForce* interatomicForce();
private:
    vector<Molecule*> m_molecules;
    vector<Atom*> m_atoms;
//    int nSimulationSteps;
    Integrator *m_integrator;

    double m_boltzmannConstant;
    double m_potentialConstant;

    mat m_boundaries;

    int m_nDimensions;

    string outFileName;
    FileFormat outFileFormat;

    mat m_cellShiftVectors;

    int pow3nDimensions;

    irowvec m_nCells;
    rowvec m_cellLengths;

    vector<MoleculeSystemCell*> m_cells;

    double m_unitLength;

    Config *m_config;

    InteratomicForce *m_interatomicForce;
};

#endif // MOLECULESYSTEM_H
