#ifndef ATOM_H
#define ATOM_H

// Local includes
class InteratomicForce;
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
    friend class InteratomicForce;
public:
    Atom(AtomType atomType);

    void setPosition(const rowvec &position);
    void setVelocity(const rowvec &velocity);
    void addForce(const rowvec &force);
    void addPotential(double potential);

    void clearForcePotentialPressure();

    void clearDisplacement();
    void addDisplacement(const rowvec &displacement);
    void addDisplacement(double displacement, uint component);

    void setCellID(int cellID);

    // Simple getters left in header for optimization
    const rowvec &position() const
    {
        return m_position;
    }
    const rowvec &velocity() const
    {
        return m_velocity;
    }
    const rowvec &force() const
    {
        return m_force;
    }
    const rowvec& displacement() const
    {
        return m_displacement;
    }
    AtomType type() const
    {
        return m_type;
    }
    double potential() const
    {
        return m_potential;
    }
    double localPressure() const {
        return m_localPressure;
    }
    int cellID() const {
        return m_cellID;
    }
    double mass()
    {
        return type().mass;
    }
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

#endif // ATOM_H
