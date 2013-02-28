#ifndef MOLECULE_H
#define MOLECULE_H

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
    Atom();

    void setPosition(const rowvec &position);
    void setVelocity(const rowvec &velocity);
    void addForce(const rowvec &force);
    void addPotential(double potential);

    void clearForcePotentialPressure();

    double mass();

    void clearDisplacement();
    void addDisplacement(const rowvec &displacement);
    void addDisplacement(double displacement, uint component);

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
    void setCellID(int cellID);
    AtomType type();
protected:
    rowvec m_position;
    rowvec m_velocity;
    rowvec m_force;
    rowvec m_displacement;
    double m_potential;
    double m_localPressure;
    AtomType m_type;

    double m_mass;
};

#endif // MOLECULE_H
