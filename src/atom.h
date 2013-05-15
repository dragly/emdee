#ifndef ATOM_H
#define ATOM_H

// Local includes
#include <src/atomtype.h>
#include <src/math/vector3.h>

// System includes
//#include <armadillo>
#include <vector>
#include <boost/serialization/serialization.hpp>

//using namespace arma;
using namespace std;

class TwoParticleForce;

/*!
 * \brief The Molecule class contains a list of atoms which are positioned
 * relatively to the molecule's position.
 */
class Atom
{
    friend class TwoParticleForce;
public:
    inline Atom();
    inline Atom(const AtomType &atomType);

    inline void clearForcePotentialPressure();

    inline void clearDisplacement();

    inline const Vector3 &position() const;
    inline const Vector3 &velocity() const;
    inline const Vector3 &force() const;
    inline const Vector3& displacement() const;
    inline const AtomType &type() const;
//    inline const int atomTypeId() const;
    inline double potential() const;
    inline double localPressure() const;
    inline int cellID() const;
    inline int id() const;
    inline double mass() const;
    inline void addLocalPressure(double pressure);
    inline void addForce(int component, double force);

//    inline void setAtomTypeId(int atomTypeId);
    inline void setPosition(const Vector3 &position);
    inline void setVelocity(const Vector3 &velocity);
    inline void addForce(const Vector3 &force);
    inline void addPotential(double potential);
    inline void setCellID(int cellID);
    inline void addDisplacement(const Vector3& displacement);
    inline void addDisplacement(double displacement, uint component);

    inline void setForce(const Vector3 &force);

    void clone(const Atom &other);
    inline void communicationClone(const Atom &other, const vector<AtomType>& moleculeSystem);

    inline bool isPositionFixed();
    inline void setPositionFixed(bool fixed);

    void setID(int id);
protected:
    Vector3 m_position;
    Vector3 m_force;
    Vector3 m_velocity;
    Vector3 m_displacement;
    double m_potential;
    double m_localPressure;

    int m_cellID;
    AtomType m_atomType;

    bool m_isPositionFixed;
    int m_id;
    int m_atomTypeId;

private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int );
};

BOOST_CLASS_IMPLEMENTATION(Atom,object_serializable)
BOOST_IS_BITWISE_SERIALIZABLE(Atom)
BOOST_IS_MPI_DATATYPE(Atom)
BOOST_CLASS_TRACKING(Atom,track_never)

// Inlined constructors
inline Atom::Atom()  :
    m_position(Vector3::createZeros()),
    m_force(Vector3::createZeros()),
    m_velocity(Vector3::createZeros()),
    m_displacement(Vector3::createZeros()),
    m_potential(0.0),
    m_localPressure(0.0),
    m_cellID(-999),
//    m_type(AtomType::argon()),
    m_isPositionFixed(false),
    m_id(-1),
    m_atomTypeId(-1)
{
}

inline Atom::Atom(const AtomType &atomType) :
    m_position(Vector3::createZeros()),
    m_force(Vector3::createZeros()),
    m_velocity(Vector3::createZeros()),
    m_displacement(Vector3::createZeros()),
    m_potential(0.0),
    m_localPressure(0.0),
    m_cellID(-999),
    m_atomType(atomType),
    m_isPositionFixed(false),
    m_id(-1),
    m_atomTypeId(atomType.id())
{
}

template<class Archive>
inline void Atom::serialize(Archive & ar, const unsigned int)
{
    ar & m_position;
    ar & m_velocity;
    ar & m_id;
    ar & m_atomTypeId;
    //        ar & m_force;
//    ar & m_displacement;
    //        ar & m_mass;
    //        ar & m_potential;
    //        ar & m_localPressure;
    //        ar & m_cellID;
    //        ar & m_type;
}

// Inlined functions
inline void Atom::addForce(int component, double force)
{
    m_force(component) += force;
}

//inline void Atom::setAtomTypeId(int atomTypeId)
//{
//    m_atomTypeId = atomTypeId;
//}

inline void Atom::addDisplacement(double displacement, uint component) {
    m_displacement(component) += displacement;
}

inline void Atom::addDisplacement(const Vector3& displacement) {
    m_displacement += displacement;
}

inline void Atom::setCellID(int cellID)
{
    m_cellID = cellID;
}
inline void Atom::addPotential(double potential) {
    m_potential += potential;
}
inline void Atom::addForce(const Vector3 &force)
{
    m_force += force;
}

inline void Atom::setForce(const Vector3 &force) {
    m_force = force;
}

inline void Atom::setID(int id) {
    m_id = id;
}

inline int Atom::id() const {
    return m_id;
}

inline void Atom::setVelocity(const Vector3 &velocity)
{
    m_velocity = velocity;
}

inline void Atom::setPosition(const Vector3 &position)
{
    m_displacement += (position - m_position);
    m_position = position;
}

inline void Atom::addLocalPressure(double pressure) {
    m_localPressure += pressure;
}

inline double Atom::mass() const
{
    return type().mass();
}

inline int Atom::cellID() const
{
    return m_cellID;
}
inline double Atom::localPressure() const
{
    return m_localPressure;
}
inline double Atom::potential() const
{
    return m_potential;
}

inline const Vector3 &Atom::position() const
{
    return m_position;
}
inline const Vector3 &Atom::velocity() const
{
    return m_velocity;
}
inline const Vector3 &Atom::force() const
{
    return m_force;
}
inline const Vector3& Atom::displacement() const
{
    return m_displacement;
}
inline const AtomType& Atom::type() const
{
    return m_atomType;
}

//inline const int Atom::atomTypeId() const
//{
//    return m_atomTypeId;
//}

inline void Atom::communicationClone(const Atom &other, const vector<AtomType>& particleTypes)
{
    this->m_id = other.m_id;
    this->m_position = other.m_position;
//    this->m_displacement = other.m_displacement;
    this->m_velocity = other.m_velocity;
    this->m_atomTypeId = other.m_atomTypeId;
    this->m_atomType = particleTypes[other.m_atomTypeId];
    //    this->m_cellID = other.m_cellID;
    //    this->m_force = other.m_force;
    //    this->m_localPressure = other.m_localPressure;
    //    this->m_potential = other.m_potential;
    //    this->m_type = other.m_type;
}

inline bool Atom::isPositionFixed() {
    return m_isPositionFixed;
}

inline void Atom::setPositionFixed(bool fixed)
{
    m_isPositionFixed = fixed;
}

inline void Atom::clearForcePotentialPressure()
{
    m_force.zeros();
    m_potential = 0;
    m_localPressure = 0;
}

inline void Atom::clearDisplacement() {
    m_displacement.zeros();
}

#endif // ATOM_H
