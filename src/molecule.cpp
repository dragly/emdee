#include "molecule.h"
#include "atom.h"

Molecule::Molecule() :
    m_position(zeros(3)),
    m_velocity(zeros(3)),
    m_force(zeros(3)),
    m_mass(1)
{
}

const vec3 &Molecule::position() const
{
    return m_position;
}


const vec3 &Molecule::velocity() const
{
    return m_velocity;
}

const vec3 &Molecule::force() const
{
    return m_force;
}

void Molecule::setPosition(const vec3 &position)
{
    m_position = position;
}


void Molecule::setVelocity(const vec3 &velocity)
{
    m_velocity = velocity;
}

void Molecule::addForce(const vec3 &force)
{
    m_force += force;
}

void Molecule::clearForces()
{
    m_force = zeros<vec>(3);
    for(Atom* atom : m_atoms) {
        atom->clearForces();
    }
}

double Molecule::mass()
{
    return m_mass;
}

void Molecule::addAtom(Atom *atom)
{
    m_atoms.push_back(atom);
}

const vector<Atom *> Molecule::atoms() {
    return m_atoms;
}
