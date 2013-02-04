#include "atom.h"
#include "molecule.h"

#include <armadillo>

using namespace arma;

Atom::Atom(Molecule *parent) :
    m_position(zeros<rowvec>(3)),
    m_velocity(zeros<rowvec>(3)),
    m_force(zeros<rowvec>(3)),
    m_parent(parent)
{
}

Atom::Atom(Molecule *parent, AtomType atomType) :
    m_position(zeros<rowvec>(3)),
    m_velocity(zeros<rowvec>(3)),
    m_force(zeros<rowvec>(3)),
    m_type(atomType),
    m_parent(parent)
{
}

void Atom::setPosition(const rowvec &position)
{
    m_position = position;
}

rowvec Atom::position()
{
    return m_position;
}

rowvec Atom::absolutePosition()
{
    if(m_parent == NULL) {
        return m_position;
    }
    return m_position + m_parent->position();
}

void Atom::setVelocity(const rowvec &velocity) {
    m_velocity = velocity;
}

const rowvec &Atom::velocity() const
{
    return m_velocity;
}

void Atom::clearForces()
{
    m_force = zeros<rowvec>(3);
}

void Atom::addForce(const rowvec &force)
{
    m_force += force;
    m_parent->addForce(force);
}

const rowvec &Atom::force() const
{
    return m_force;
}

rowvec Atom::absoluteVelocity() {
    if(m_parent == NULL) {
        return m_velocity;
    }
    return m_velocity + m_parent->velocity();
}

double Atom::mass()
{
    return type().mass;
}

AtomType Atom::type()
{
    return m_type;
}
