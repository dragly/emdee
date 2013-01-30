#ifndef MOLECULESYSTEM_H
#define MOLECULESYSTEM_H

// Forward declarations
class Atom;
class Molecule;
class Integrator;

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

    void load(string fileName);
    bool save(string fileName);
    bool saveXyz(string fileName);
    void addMolecules(vector<Molecule *> molecule);
    const vector<Molecule*> &molecules() const;
    void updateForces();
    void simulate();

    void setBoundaries(double min, double max);
    void setBoundaries(double xMin, double xMax, double yMin, double yMax, double zMin, double zMax);
    void setBoundaries(mat::fixed<2, 3> boundaries);

private:
    vector<Molecule*> m_molecules;
    vector<Atom*> m_atoms;
    int nSimulationSteps;
    Integrator *integrator;

    double boltzmannConstant;
    double potentialConstant;

    mat::fixed<2,3> m_boundaries;

    int nDimensions;
};

#endif // MOLECULESYSTEM_H
