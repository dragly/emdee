#include <src/atom.h>
#include <src/molecule.h>

#include <armadillo>

using namespace arma;

Atom::Atom(Molecule *parent) :
    m_relativePosition(zeros<rowvec>(3)),
    m_relativeVelocity(zeros<rowvec>(3)),
    m_force(zeros<rowvec>(3)),
    m_potential(0),
    m_parent(parent),
    m_localPressure(0),
    m_cellID(-999)
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

void Atom::setRelativeVelocity(const rowvec &velocity) {
    refreshAbsolutePositionAndVelocity();
    m_relativeVelocity = velocity;
}

const rowvec &Atom::relativeVelocity() const
{
    return m_relativeVelocity;
}

void Atom::clearForcePotentialPressure()
{
    m_force.zeros();
    m_potential = 0;
    m_localPressure = 0;
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

double Atom::mass()
{
    return type().mass;
}

AtomType Atom::type()
{
    return m_type;
}

void Atom::setCellID(int cellID)
{
    m_cellID = cellID;
}

const rowvec& Atom::displacement() const {
    return m_parent->displacement();
}
