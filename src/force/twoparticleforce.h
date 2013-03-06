#ifndef TWOPARTICLEFORCE_H
#define TWOPARTICLEFORCE_H

// Local headers
#include <src/math/vector3.h>

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
    virtual void calculateAndApplyForce(Atom* atom1, Atom* atom2, const Vector3 &atom2Offset) = 0;

    void setNewtonsThirdLawEnabled(bool enable);
    bool isNewtonsThirdLawEnabled();
protected:
    bool m_isNewtonsThirdLawEnabled;
};

inline void TwoParticleForce::setNewtonsThirdLawEnabled(bool enable) {
    m_isNewtonsThirdLawEnabled = enable;
}

inline bool TwoParticleForce::isNewtonsThirdLawEnabled() {
    return m_isNewtonsThirdLawEnabled;
}

#endif // TWOPARTICLEFORCE_H
