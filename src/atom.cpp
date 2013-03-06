#include "atom.h"

//void Atom::


//void Atom::

void Atom::clearForcePotentialPressure()
{
    m_force.zeros();
    m_potential = 0;
    m_localPressure = 0;
}

void Atom::clearDisplacement() {
    m_displacement.zeros();
}

void Atom::clone(const Atom &other)
{
    this->m_position = other.m_position;
    this->m_cellID = other.m_cellID;
    this->m_displacement = other.m_displacement;
    this->m_force = other.m_force;
    this->m_localPressure = other.m_localPressure;
    this->m_potential = other.m_potential;
    this->m_type = other.m_type;
    this->m_velocity = other.m_velocity;
}

//void Atom::

//void Atom::

//void Atom::

//void Atom::

//void Atom::

//void Atom::

//void Atom::

//template<class Archive>
//void Atom::serialize(Archive & ar, const unsigned int version)
//{
//    ar & m_position;
//    ar & m_velocity;
//    ar & m_force;
//    ar & m_displacement;
//    ar & m_mass;
//    ar & m_potential;
//    ar & m_localPressure;
//    ar & m_cellID;
//    ar & m_type;
//}
