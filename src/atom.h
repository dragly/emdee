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

    void setPosition(const vec &position);
    vec position();
    void setVelocity(const vec &velocity);
    const vec &velocity() const;
    void addForce(const vec &force);
    const vec &force() const;
    AtomType type();

    vec absolutePosition();
    vec absoluteVelocity();

    double mass();
    void clearForces();
private:
    vec m_position;
    vec m_velocity;
    vec m_force;
    AtomType m_type;

    Molecule* m_parent;
};

#endif // ATOM_H
