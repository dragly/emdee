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

void Atom::refreshAbsolutePositionAndVelocity() {
    m_absolutePosition = m_position + m_parent->position();
    m_absoluteVelocity = m_velocity + m_parent->velocity();
}

void Atom::setPosition(const rowvec &position)
{
    refreshAbsolutePositionAndVelocity();
    m_position = position;
}

const rowvec& Atom::position() const
{
    return m_position;
}

const rowvec& Atom::absolutePosition() const
{
//    if(m_parent == NULL) {
//        return m_position;
//    }
    return m_absolutePosition;
}

void Atom::setVelocity(const rowvec &velocity) {
    refreshAbsolutePositionAndVelocity();
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

const rowvec& Atom::absoluteVelocity() const
{
//    if(m_parent == NULL) {
//        return m_velocity;
//    }
    return m_absoluteVelocity;
}

double Atom::mass()
{
    return type().mass;
}

AtomType Atom::type()
{
    return m_type;
}
