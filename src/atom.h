#ifndef ATOM_H
#define ATOM_H

// Local includes
class TwoParticleForce;
#include <src/atomtype.h>

// System includes
#include <armadillo>
#include <vector>

using namespace arma;
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
    const rowvec &position() const;
    const rowvec &velocity() const;
    const rowvec &force() const;
    const rowvec& displacement() const;
    AtomType type() const;
    double potential() const;
    double localPressure() const;
    int cellID() const;
    double mass() const;
    void addLocalPressure(double pressure);
    void addForce(int component, double force);

    void setPosition(const rowvec &position);
    void setVelocity(const rowvec &velocity);
    void addForce(const rowvec &force);
    void addPotential(double potential);
    void setCellID(int cellID);
    void addDisplacement(const rowvec& displacement);
    void addDisplacement(double displacement, uint component);

protected:
    rowvec m_position;
    rowvec m_velocity;
    rowvec m_force;
    rowvec m_displacement;
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

inline void Atom::addDisplacement(const rowvec& displacement) {
    m_displacement += displacement;
}

inline void Atom::setCellID(int cellID)
{
    m_cellID = cellID;
}
inline void Atom::addPotential(double potential) {
    m_potential += potential;
}
inline void Atom::addForce(const rowvec &force)
{
    m_force += force;
}
inline void Atom::setVelocity(const rowvec &velocity)
{
    m_velocity = velocity;
}
inline void Atom::setPosition(const rowvec &position)
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

inline const rowvec &Atom::position() const
{
    return m_position;
}
inline const rowvec &Atom::velocity() const
{
    return m_velocity;
}
inline const rowvec &Atom::force() const
{
    return m_force;
}
inline const rowvec& Atom::displacement() const
{
    return m_displacement;
}
inline AtomType Atom::type() const
{
    return m_type;
}

#endif // ATOM_H
