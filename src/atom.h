#ifndef ATOM_H
#define ATOM_H

class Molecule;
#include "atomtype.h"

#include <armadillo>

using namespace arma;
using namespace std;

class Atom
{
public:
    Atom(Molecule* parent);
    Atom(Molecule* parent, AtomType atomType);

    void setPosition(const vec3 &position);
    vec3 position();
    void setVelocity(const vec3 &velocity);
    const vec3 &velocity() const;
    void addForce(const vec3 &force);
    const vec3 &force() const;
    AtomType type();

    vec3 absolutePosition();
    vec3 absoluteVelocity();

    double mass();
    void clearForces();
private:
    vec3 m_position;
    vec3 m_velocity;
    vec3 m_force;
    AtomType m_type;

    Molecule* m_parent;
};

#endif // ATOM_H
