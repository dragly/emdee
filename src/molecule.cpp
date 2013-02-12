#include <src/molecule.h>
#include <src/atom.h>

Molecule::Molecule() :
    m_position(zeros<rowvec>(3)),
    m_velocity(zeros<rowvec>(3)),
    m_force(zeros<rowvec>(3)),
    m_mass(0.0)
{
}

const rowvec &Molecule::position() const
{
    return m_position;
}


const rowvec &Molecule::velocity() const
{
    return m_velocity;
}

const rowvec &Molecule::force() const
{
    return m_force;
}

void Molecule::setPosition(const rowvec &position)
{
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
