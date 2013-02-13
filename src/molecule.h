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

    const rowvec &position() const;
    const rowvec &velocity() const;
    const rowvec &force() const;

    void setPosition(const rowvec &position);
    void setVelocity(const rowvec &velocity);
    void addForce(const rowvec &force);

    void clearForces();

    double mass();

    void addAtom(Atom* atom);

    const vector<Atom *> atoms();
protected:
    rowvec m_position;
    rowvec m_velocity;
    rowvec m_force;

    double m_mass;

    vector<Atom*> m_atoms;
};

#endif // MOLECULE_H
