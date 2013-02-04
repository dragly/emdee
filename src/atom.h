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

    void setPosition(const rowvec &position);
    rowvec position();
    void setVelocity(const rowvec &velocity);
    const rowvec &velocity() const;
    void addForce(const rowvec &force);
    const rowvec &force() const;
    AtomType type();

    rowvec absolutePosition();
    rowvec absoluteVelocity();

    double mass();
    void clearForces();
private:
    rowvec m_position;
    rowvec m_velocity;
    rowvec m_force;
    AtomType m_type;

    Molecule* m_parent;
};

#endif // ATOM_H
