#ifndef ATOM_H
#define ATOM_H

class Molecule;
#include <src/atomtype.h>

#include <armadillo>

using namespace arma;
using namespace std;

class Atom
{
public:
    Atom(Molecule* parent);
    Atom(Molecule* parent, AtomType atomType);

    void setRelativePosition(const rowvec &position);
    const rowvec& relativePosition() const;
    void setRelativeVelocity(const rowvec &velocity);
    const rowvec &relativeVelocity() const;
    void addForce(const rowvec &force);
    void addPotential(double potential);
    const rowvec &force() const;
    AtomType type();

    const rowvec& position() const;
    const rowvec& velocity() const;

    double mass();
    void clearForceAndPotential();
    void refreshAbsolutePositionAndVelocity();

    double potential();
protected:
    rowvec m_relativePosition;
    rowvec m_relativeVelocity;
    rowvec m_force;
    rowvec m_position;
    rowvec m_velocity;
    AtomType m_type;
    double m_potential;

    Molecule* m_parent;
};

#endif // ATOM_H
