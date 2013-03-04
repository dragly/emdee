#include "atom.h"

Atom::Atom()  :
    m_position(Vector3::createZeros()),
    m_velocity(Vector3::createZeros()),
    m_force(Vector3::createZeros()),
    m_displacement(Vector3::createZeros()),
    m_mass(0.0),
    m_potential(0.0),
    m_localPressure(0.0),
    m_cellID(-999),
    m_type(AtomType::argon())
{

}

Atom::Atom(AtomType atomType) :
    m_position(Vector3::createZeros()),
    m_velocity(Vector3::createZeros()),
    m_force(Vector3::createZeros()),
    m_displacement(Vector3::createZeros()),
    m_mass(0.0),
    m_potential(0.0),
    m_localPressure(0.0),
    m_cellID(-999),
    m_type(atomType)
{
}

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
