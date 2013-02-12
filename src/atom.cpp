#include <src/atom.h>
#include <src/molecule.h>

#include <armadillo>

using namespace arma;

Atom::Atom(Molecule *parent) :
    m_relativePosition(zeros<rowvec>(3)),
    m_relativeVelocity(zeros<rowvec>(3)),
    m_force(zeros<rowvec>(3)),
    m_potential(0),
    m_parent(parent)
{
}

Atom::Atom(Molecule *parent, AtomType atomType) :
    m_relativePosition(zeros<rowvec>(3)),
    m_relativeVelocity(zeros<rowvec>(3)),
    m_force(zeros<rowvec>(3)),
    m_type(atomType),
    m_parent(parent)
{
}

void Atom::refreshAbsolutePositionAndVelocity() {
    m_position = m_relativePosition + m_parent->position();
    m_velocity = m_relativeVelocity + m_parent->velocity();
}

double Atom::potential()
{
    return m_potential;
}

void Atom::addPotential(double potential) {
    m_potential += potential;
}

void Atom::setRelativePosition(const rowvec &position)
{
    refreshAbsolutePositionAndVelocity();
    m_relativePosition = position;
}

const rowvec& Atom::relativePosition() const
{
    return m_relativePosition;
}

const rowvec& Atom::position() const
{
//    if(m_parent == NULL) {
//        return m_position;
//    }
    return m_position;
}

void Atom::setRelativeVelocity(const rowvec &velocity) {
    refreshAbsolutePositionAndVelocity();
    m_relativeVelocity = velocity;
}

const rowvec &Atom::relativeVelocity() const
{
    return m_relativeVelocity;
}

void Atom::clearForceAndPotential()
{
    m_force.zeros();
    m_potential = 0;
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

const rowvec& Atom::velocity() const
{
//    if(m_parent == NULL) {
//        return m_velocity;
//    }
    return m_velocity;
}

double Atom::mass()
{
    return type().mass;
}

AtomType Atom::type()
{
    return m_type;
}
