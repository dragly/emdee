#include <src/atom.h>

Atom::Atom(AtomType atomType) :
    m_position(zeros<rowvec>(3)),
    m_velocity(zeros<rowvec>(3)),
    m_force(zeros<rowvec>(3)),
    m_displacement(zeros<rowvec>(3)),
    m_mass(0.0),
    m_potential(0.0),
    m_localPressure(0.0),
    m_cellID(-999),
    m_type(atomType)
{
}

void Atom::setPosition(const rowvec &position)
{
    m_displacement += (position - m_position);
    m_position = position;
}


void Atom::setVelocity(const rowvec &velocity)
{
    m_velocity = velocity;
}

void Atom::addForce(const rowvec &force)
{
    m_force += force;
}

void Atom::clearForcePotentialPressure()
{
    m_force = zeros<rowvec>(3);
    m_potential = 0;
    m_localPressure = 0;
}

void Atom::clearDisplacement() {
    m_displacement.zeros();
}

void Atom::addDisplacement(const rowvec& displacement) {
    m_displacement += displacement;
}

void Atom::addDisplacement(double displacement, uint component) {
    m_displacement(component) += displacement;
}

void Atom::addPotential(double potential) {
    m_potential += potential;
}


void Atom::setCellID(int cellID)
{
    m_cellID = cellID;
}
