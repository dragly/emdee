#ifndef ATOM_H
#define ATOM_H

// Local includes
class TwoParticleForce;
#include <src/atomtype.h>
#include <src/math/vector3.h>

// System includes
//#include <armadillo>
#include <vector>

//using namespace arma;
using namespace std;

/*!
 * \brief The Molecule class contains a list of atoms which are positioned
 * relatively to the molecule's position.
 */
class Atom
{
    friend class TwoParticleForce;
public:
    Atom(AtomType atomType);

    void clearForcePotentialPressure();

    void clearDisplacement();

    // Simple getters left in header for optimization
    const Vector3 &position() const;
    const Vector3 &velocity() const;
    const Vector3 &force() const;
    const Vector3& displacement() const;
    AtomType type() const;
    double potential() const;
    double localPressure() const;
    int cellID() const;
    double mass() const;
    void addLocalPressure(double pressure);
    void addForce(int component, double force);

    void setPosition(const Vector3 &position);
    void setVelocity(const Vector3 &velocity);
    void addForce(const Vector3 &force);
    void addPotential(double potential);
    void setCellID(int cellID);
    void addDisplacement(const Vector3& displacement);
    void addDisplacement(double displacement, uint component);

protected:
    Vector3 m_position;
    Vector3 m_velocity;
    Vector3 m_force;
    Vector3 m_displacement;
    double m_mass;
    double m_potential;
    double m_localPressure;

    int m_cellID;
    AtomType m_type;
};

// Inlined functions
inline void Atom::addForce(int component, double force)
{
   m_force(component) += force;
}

inline void Atom::addDisplacement(double displacement, uint component) {
    m_displacement(component) += displacement;
}

inline void Atom::addDisplacement(const Vector3& displacement) {
    m_displacement += displacement;
}

inline void Atom::setCellID(int cellID)
{
    m_cellID = cellID;
}
inline void Atom::addPotential(double potential) {
    m_potential += potential;
}
inline void Atom::addForce(const Vector3 &force)
{
    m_force += force;
}
inline void Atom::setVelocity(const Vector3 &velocity)
{
    m_velocity = velocity;
}
inline void Atom::setPosition(const Vector3 &position)
{
    m_displacement += (position - m_position);
    m_position = position;
}
inline void Atom::addLocalPressure(double pressure) {
    m_localPressure += pressure;
}
inline double Atom::mass() const
{
    return type().mass;
}
inline int Atom::cellID() const
{
    return m_cellID;
}
inline double Atom::localPressure() const
{
    return m_localPressure;
}
inline double Atom::potential() const
{
    return m_potential;
}

inline const Vector3 &Atom::position() const
{
    return m_position;
}
inline const Vector3 &Atom::velocity() const
{
    return m_velocity;
}
inline const Vector3 &Atom::force() const
{
    return m_force;
}
inline const Vector3& Atom::displacement() const
{
    return m_displacement;
}
inline AtomType Atom::type() const
{
    return m_type;
}

#endif // ATOM_H
