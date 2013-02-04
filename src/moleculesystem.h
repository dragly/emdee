#ifndef MOLECULESYSTEM_H
#define MOLECULESYSTEM_H

// Forward declarations
class Atom;
class Molecule;
class Integrator;
class MoleculeSystemCell;

// System includes
#include <iostream>
#include <vector>
#include <armadillo>

using namespace std;
using namespace arma;

class MoleculeSystem
{
public:
    MoleculeSystem();

    enum FileFormat {
        XyzFormat,
        HDF5Format
    };

    void load(string fileName);
    bool save(int step);
    bool saveXyz(int step);
    void addMolecules(const vector<Molecule *> &molecule);
    const vector<Molecule*> &molecules() const;
    void updateForces();
    void simulate();

    void setBoundaries(double min, double max);
    void setBoundaries(double xMin, double xMax, double yMin, double yMax, double zMin, double zMax);
    void setBoundaries(mat boundaries);

    bool saveHDF5(string filename);
    void setupCells(double minCutLength);

    double potentialConstant();

    const mat& cellShiftVectors();

    void refreshCellContents();
private:
    vector<Molecule*> m_molecules;
    vector<Atom*> m_atoms;
    int nSimulationSteps;
    Integrator *integrator;

    double boltzmannConstant;
    double m_potentialConstant;

    mat m_boundaries;

    int nDimensions;

    string outFileName;
    FileFormat outFileFormat;

    mat m_cellShiftVectors;

    int pow3nDimensions;

    irowvec m_nCells;
    rowvec m_cellLengths;

    vector<MoleculeSystemCell*> m_cells;
};

#endif // MOLECULESYSTEM_H
