#ifndef THREEPARTICLEFORCE_H
#define THREEPARTICLEFORCE_H

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
class ThreeParticleForce
{
public:
    ThreeParticleForce();

    virtual void calculateAndApplyForce(Atom *atom1, Atom *atom2, Atom *atom3) = 0;
    virtual void calculateAndApplyForce(Atom* atom1, Atom* atom2, Atom *atom3, const Vector3 &atom2Offset, const Vector3 &atom3Offset) = 0;

    void setNewtonsThirdLawEnabled(bool enable);
    bool isNewtonsThirdLawEnabled();
    void setCalculatePressureEnabled(bool enable);
    bool isCalculatePressureEnabled();
    void setCalculatePotentialEnabled(bool enable);
    bool isCalculatePotentialEnabled();
protected:
    bool m_isNewtonsThirdLawEnabled;
    bool m_isCalculatePressureEnabled;
    bool m_isCalculatePotentialEnabled;
    Vector3 m_zeroVector;
};

inline void ThreeParticleForce::setNewtonsThirdLawEnabled(bool enable) {
    m_isNewtonsThirdLawEnabled = enable;
}

inline bool ThreeParticleForce::isNewtonsThirdLawEnabled() {
    return m_isNewtonsThirdLawEnabled;
}

inline void ThreeParticleForce::setCalculatePressureEnabled(bool enable) {
    m_isCalculatePressureEnabled = enable;
}

inline bool ThreeParticleForce::isCalculatePressureEnabled() {
    return m_isCalculatePressureEnabled;
}

inline void ThreeParticleForce::setCalculatePotentialEnabled(bool enable) {
    m_isCalculatePotentialEnabled = enable;
}

inline bool ThreeParticleForce::isCalculatePotentialEnabled() {
    return m_isCalculatePotentialEnabled;
}

#endif // THREEPARTICLEFORCE_H
