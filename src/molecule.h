#ifndef MOLECULE_H
#define MOLECULE_H

// Local includes
class Atom;

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
public:
    Molecule();

    const vec &position() const;
    const vec &velocity() const;
    const vec &force() const;

    void setPosition(const vec &position);
    void setVelocity(const vec &velocity);
    void addForce(const vec &force);

    void clearForces();

    double mass();

    void addAtom(Atom* atom);

    const vector<Atom *> atoms();
private:
    vec m_position;
    vec m_velocity;
    vec m_force;

    double m_mass;

    vector<Atom*> m_atoms;
};

#endif // MOLECULE_H
