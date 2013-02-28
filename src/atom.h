#ifndef ATOM_H
#define ATOM_H

class Atom;
class InteratomicForce;
#include <src/atomtype.h>

#include <armadillo>

using namespace arma;
using namespace std;

class Atom_old
{
    friend class InteratomicForce;
public:
    Atom_old(Atom* parent);
    Atom_old(Atom* parent, AtomType atomType);

    void setRelativePosition(const rowvec &position);
    const rowvec& relativePosition() const;
    void setRelativeVelocity(const rowvec &velocity);
    const rowvec &relativeVelocity() const;
    void addForce(const rowvec &force);
    void addPotential(double potential);
    const rowvec &force() const;
    AtomType type();
    void setCellID(int cellID);

    double mass();
    void clearForcePotentialPressure();
    void refreshAbsolutePositionAndVelocity();

    double potential();

    // Left in header for optimization
    const rowvec& velocity() const
    {
        return m_velocity;
    }
    const rowvec& position() const
    {
        return m_position;
    }
    int cellID() {
        return m_cellID;
    }
    double localPressure() const {
        return m_localPressure;
    }

    const rowvec& displacement() const;
protected:
    rowvec m_relativePosition;
    rowvec m_relativeVelocity;
    rowvec m_position;
    rowvec m_force;
    rowvec m_velocity;
    AtomType m_type;
    double m_potential;
    double m_localPressure;

    Atom* m_parent;

    int m_cellID;
};

#endif // ATOM_H
