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

void Atom::setPosition(const vec3 &position)
{
    m_position = position;
}

vec3 Atom::position()
{
    return m_position;
}

vec3 Atom::absolutePosition()
{
    if(m_parent == NULL) {
        return m_position;
    }
    return m_position + m_parent->position();
}

void Atom::setVelocity(const vec3 &velocity) {
    m_velocity = velocity;
}

const vec3 &Atom::velocity() const
{
    return m_velocity;
}

void Atom::clearForces()
{
    m_force = zeros<vec>(3);
}

void Atom::addForce(const vec3 &force)
{
    m_force += force;
    m_parent->addForce(force);
}

const vec3 &Atom::force() const
{
    return m_force;
}

vec3 Atom::absoluteVelocity() {
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
