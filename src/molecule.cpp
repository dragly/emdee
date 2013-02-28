#include <src/molecule.h>
#include <src/atom.h>

Atom::Atom() :
    m_position(zeros<rowvec>(3)),
    m_velocity(zeros<rowvec>(3)),
    m_force(zeros<rowvec>(3)),
    m_displacement(zeros<rowvec>(3)),
    m_mass(0.0)
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
}

double Atom::mass()
{
    return m_mass;
}

//void Atom::addAtom(Atom_old *atom)
//{
//    m_atoms.push_back(atom);
//    m_mass += atom->mass();
//}

//const vector<Atom_old *> Atom::atoms() {
//    return m_atoms;
//}

void Atom::clearDisplacement() {
    m_displacement.zeros();
}

void Atom::addDisplacement(const rowvec& displacement) {
    m_displacement += displacement;
}

void Atom::addDisplacement(double displacement, uint component) {
    m_displacement(component) += displacement;
}
