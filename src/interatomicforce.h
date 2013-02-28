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

/*!
 * \brief The InteratomicForce class calculates forces between atoms.
 */
class InteratomicForce
{
public:
    InteratomicForce();

    void calculateAndApplyForce(Atom *atom1, Atom *atom2);
    void calculateAndApplyForce(Atom* atom1, Atom* atom2, const rowvec &atom2Offset);

    void setPotentialConstant(double potentialConstant);
    void setEnergyConstant(double energyConstant);

protected:
    rowvec zeroVector;

    double m_potentialConstant;
    double m_potentialConstantSquared;
    double m_energyConstant;
    double m_energyConstant4;
    double m_energyConstant24;
};

#endif // INTERATOMICFORCE_H
