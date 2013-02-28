#include <src/atom.h>
#include <src/molecule.h>

#include <armadillo>

using namespace arma;

Atom_old::Atom_old(Atom *parent) :
    m_relativePosition(zeros<rowvec>(3)),
    m_relativeVelocity(zeros<rowvec>(3)),
    m_force(zeros<rowvec>(3)),
    m_potential(0),
    m_localPressure(0),
    m_parent(parent),
    m_cellID(-999)
{
}

Atom_old::Atom_old(Atom *parent, AtomType atomType) :
    m_relativePosition(zeros<rowvec>(3)),
    m_relativeVelocity(zeros<rowvec>(3)),
    m_force(zeros<rowvec>(3)),
    m_type(atomType),
    m_parent(parent)
{
}

void Atom_old::refreshAbsolutePositionAndVelocity() {
    m_position = m_relativePosition + m_parent->position();
    m_velocity = m_relativeVelocity + m_parent->velocity();
}

double Atom_old::potential()
{
    return m_potential;
}

void Atom_old::addPotential(double potential) {
    m_potential += potential;
}

void Atom_old::setRelativePosition(const rowvec &position)
{
    refreshAbsolutePositionAndVelocity();
    m_relativePosition = position;
}

const rowvec& Atom_old::relativePosition() const
{
    return m_relativePosition;
}

void Atom_old::setRelativeVelocity(const rowvec &velocity) {
    refreshAbsolutePositionAndVelocity();
    m_relativeVelocity = velocity;
}

const rowvec &Atom_old::relativeVelocity() const
{
    return m_relativeVelocity;
}

void Atom_old::clearForcePotentialPressure()
{
    m_force.zeros();
    m_potential = 0;
    m_localPressure = 0;
}

void Atom_old::addForce(const rowvec &force)
{
    m_force += force;
    m_parent->addForce(force);
}

const rowvec &Atom_old::force() const
{
    return m_force;
}

double Atom_old::mass()
{
    return type().mass;
}

AtomType Atom_old::type()
{
    return m_type;
}

void Atom_old::setCellID(int cellID)
{
    m_cellID = cellID;
}

const rowvec& Atom_old::displacement() const {
    return m_parent->displacement();
}
