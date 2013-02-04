#include "atom.h"
#include "molecule.h"

#include <armadillo>

using namespace arma;

Atom::Atom(Molecule *parent) :
    m_position(zeros(3)),
    m_velocity(zeros(3)),
    m_force(zeros(3)),
    m_parent(parent)
{
}

Atom::Atom(Molecule *parent, AtomType atomType) :
    m_position(zeros(3)),
    m_velocity(zeros(3)),
    m_force(zeros(3)),
    m_type(atomType),
    m_parent(parent)
{
}

void Atom::setPosition(const vec &position)
{
    m_position = position;
}

vec Atom::position()
{
    return m_position;
}

vec Atom::absolutePosition()
{
    if(m_parent == NULL) {
        return m_position;
    }
    return m_position + m_parent->position();
}

void Atom::setVelocity(const vec &velocity) {
    m_velocity = velocity;
}

const vec &Atom::velocity() const
{
    return m_velocity;
}

void Atom::clearForces()
{
    m_force = zeros<vec>(3);
}

void Atom::addForce(const vec &force)
{
    m_force += force;
    m_parent->addForce(force);
}

const vec &Atom::force() const
{
    return m_force;
}

vec Atom::absoluteVelocity() {
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
