#ifndef INTERATOMICFORCE_H
#define INTERATOMICFORCE_H

// Local headers

// Forward declarations
class Atom;
class MoleculeSystem;

// System headers
#include <armadillo>
#include <libconfig.h++>

using namespace arma;
using namespace libconfig;


class InteratomicForce
{
public:
    InteratomicForce(MoleculeSystem *parent);

    void loadConfiguration(Config* config);

    const rowvec force(Atom *atom1, Atom *atom2);
    const rowvec force(Atom* atom1, Atom* atom2, const rowvec &atom2Offset);
private:
    rowvec tmpForce;
    rowvec zeroVector;

    double m_potentialConstant;

    MoleculeSystem* moleculeSystem;
    rowvec rVec;
};

#endif // INTERATOMICFORCE_H
