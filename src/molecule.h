#ifndef MOLECULE_H
#define MOLECULE_H

// Local includes
class Atom;
class InteratomicForce;

// System includes
#include <armadillo>
#include <vector>

using namespace arma;
using namespace std;

/*!
 * \brief The Molecule class contains a list of atoms which are positioned
 * relatively to the molecule's position.
 */
class Molecule
{
    friend class InteratomicForce;
public:
    Molecule();

    void setPosition(const rowvec &position);
    void setVelocity(const rowvec &velocity);
    void addForce(const rowvec &force);

    void clearForces();

    double mass();

    void addAtom(Atom* atom);

    const vector<Atom *> atoms();
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
protected:
    rowvec m_position;
    rowvec m_velocity;
    rowvec m_force;
    rowvec m_displacement;

    double m_mass;

    vector<Atom*> m_atoms;
};

#endif // MOLECULE_H
