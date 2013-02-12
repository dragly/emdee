#ifndef GENERATOR_H
#define GENERATOR_H

class Molecule;
class Atom;

#include <src/atomtype.h>

#include <vector>
#include <armadillo>
#include <libconfig.h++>

using namespace std;
using namespace arma;
using namespace libconfig;

class Generator
{
public:
    Generator();
//    void loadConfiguration(Config* config);
    vector<Molecule*> generateFcc(double b, int nCells, AtomType atomType);
    void boltzmannDistributeVelocities(double temperature, vector<Molecule *> molecules);
    void uniformDistributeVelocities(double maxVelocity, vector<Molecule *> molecules);

//    void setUnitLength(double unitLength);

    const mat& lastBoundaries() const;

    void setNCells(int nCells);
    void setB(double b);
    void setBoltzmannConstant(double boltzmannConstant);
protected:
    double m_unitLength;
    double m_boltzmannConstant;

    mat m_lastBoundaries;

    int m_nDimensions;
};

#endif // GENERATOR_H
