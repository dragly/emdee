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

    const vec3 &position() const;
    const vec3 &velocity() const;
    const vec3 &force() const;

    void setPosition(const vec3 &position);
    void setVelocity(const vec3 &velocity);
    void addForce(const vec3 &force);

    void clearForces();

    double mass();

    void addAtom(Atom* atom);

    const vector<Atom *> atoms();
private:
    vec3 m_position;
    vec3 m_velocity;
    vec3 m_force;

    double m_mass;

    vector<Atom*> m_atoms;
};

#endif // MOLECULE_H
