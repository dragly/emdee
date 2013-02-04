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
    const rowvec& position() const;
    void setVelocity(const rowvec &velocity);
    const rowvec &velocity() const;
    void addForce(const rowvec &force);
    const rowvec &force() const;
    AtomType type();

    const rowvec& absolutePosition() const;
    const rowvec& absoluteVelocity() const;

    double mass();
    void clearForces();
    void refreshAbsolutePositionAndVelocity();
private:
    rowvec m_position;
    rowvec m_velocity;
    rowvec m_force;
    rowvec m_absolutePosition;
    rowvec m_absoluteVelocity;
    AtomType m_type;

    Molecule* m_parent;
};

#endif // ATOM_H
