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
    InteratomicForce();

//    void loadConfiguration(Config* config);

    void calculateAndApplyForce(Atom *atom1, Atom *atom2);
    void calculateAndApplyForce(Atom* atom1, Atom* atom2, const rowvec &atom2Offset);

//    const rowvec& force();
//    double potential();

    void setPotentialConstant(double potentialConstant);
    void setEnergyConstant(double energyConstant);

protected:
    rowvec tmpForce;
    double tmpPotential;
    rowvec zeroVector;

    double m_potentialConstant;
    double m_energyConstant;
};

#endif // INTERATOMICFORCE_H
