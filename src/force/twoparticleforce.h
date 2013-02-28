#ifndef TWOPARTICLEFORCE_H
#define TWOPARTICLEFORCE_H

// Local headers

// Forward declarations
class Atom;

// System headers
#include <armadillo>
#include <libconfig.h++>

using namespace arma;
using namespace libconfig;

/*!
 * \brief The InteratomicForce class calculates forces between atoms.
 */
class TwoParticleForce
{
public:
    TwoParticleForce();

    virtual void calculateAndApplyForce(Atom *atom1, Atom *atom2) = 0;
    virtual void calculateAndApplyForce(Atom* atom1, Atom* atom2, const rowvec &atom2Offset) = 0;

protected:
};

#endif // TWOPARTICLEFORCE_H
