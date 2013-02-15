#include <src/molecule.h>
#include <src/atom.h>

Molecule::Molecule() :
    m_position(zeros<rowvec>(3)),
    m_velocity(zeros<rowvec>(3)),
    m_force(zeros<rowvec>(3)),
    m_displacement(zeros<rowvec>(3)),
    m_mass(0.0)
{
}

void Molecule::setPosition(const rowvec &position)
{
    m_displacement += (position - m_position);
    m_position = position;
    for(Atom* atom : m_atoms) {
        atom->refreshAbsolutePositionAndVelocity();
    }
}


void Molecule::setVelocity(const rowvec &velocity)
{
    m_velocity = velocity;
    for(Atom* atom : m_atoms) {
        atom->refreshAbsolutePositionAndVelocity();
    }
}

void Molecule::addForce(const rowvec &force)
{
    m_force += force;
}

void Molecule::clearForces()
{
    m_force = zeros<rowvec>(3);
    for(Atom* atom : m_atoms) {
        atom->clearForceAndPotential();
    }
}

double Molecule::mass()
{
    return m_mass;
}

void Molecule::addAtom(Atom *atom)
{
    m_atoms.push_back(atom);
    m_mass += atom->mass();
}

const vector<Atom *> Molecule::atoms() {
    return m_atoms;
}

void Molecule::clearDisplacement() {
    m_displacement.zeros();
}

void Molecule::addDisplacement(const rowvec& displacement) {
    m_displacement += displacement;
}

void Molecule::addDisplacement(double displacement, uint component) {
    m_displacement(component) += displacement;
}
