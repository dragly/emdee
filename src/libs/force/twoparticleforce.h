#ifndef TWOPARTICLEFORCE_H
#define TWOPARTICLEFORCE_H

// Local headers
#include <math/vector3.h>

// Forward declarations
class Atom;

// System headers
#include <armadillo>
//#include <libconfig.h++>

using namespace arma;
//using namespace libconfig;

/*!
 * \brief The InteratomicForce class calculates forces between atoms.
 */
class TwoParticleForce
{
public:
    TwoParticleForce();

    virtual void calculateAndApplyForce(Atom *atom1, Atom *atom2);
    virtual void calculateAndApplyForce(Atom* atom1, Atom* atom2, const Vector3 &atom2Offset) = 0;

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

inline void TwoParticleForce::setNewtonsThirdLawEnabled(bool enable) {
    m_isNewtonsThirdLawEnabled = enable;
}

inline bool TwoParticleForce::isNewtonsThirdLawEnabled() {
    return m_isNewtonsThirdLawEnabled;
}

inline void TwoParticleForce::setCalculatePressureEnabled(bool enable) {
    m_isCalculatePressureEnabled = enable;
}

inline bool TwoParticleForce::isCalculatePressureEnabled() {
    return m_isCalculatePressureEnabled;
}

inline void TwoParticleForce::setCalculatePotentialEnabled(bool enable) {
    m_isCalculatePotentialEnabled = enable;
}

inline bool TwoParticleForce::isCalculatePotentialEnabled() {
    return m_isCalculatePotentialEnabled;
}

#endif // TWOPARTICLEFORCE_H
