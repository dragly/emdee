#ifndef ATOMICTYPE_H
#define ATOMICTYPE_H

#include <string>

#ifdef USE_MPI
#include <boost/mpi.hpp>
#include <boost/serialization/string.hpp>
#endif

using namespace std;

class AtomType
{
public:
    enum AtomTypeEnum {
        HYDROGEN,
        HELIUM,
        OXYGEN,
        ARGON
    };

    AtomType();
    AtomType(int index);
    AtomType(AtomTypeEnum atomTypeEnum, int index);
    static AtomType argon();

    void setName(string name);
    void setAbbreviation(string abbreviation);
    void setNumber(int number);
    void setMass(double mass);

    int number() const;

    string name() const;

    string abbreviation() const;

    double effectiveCharge() const;

    double electronicPolarizability() const;
    void setEffectiveCharge(double arg);

    void setElectronicPolarizability(double arg);

    double mass() const;

    int index() const;

    bool operator ==(const AtomType& other) const;
protected:
    string m_name;
    string m_abbreviation;
    int m_number;
    double m_mass;
    double m_effectiveCharge;
    double m_electronicPolarizability;
    int m_index;

//private:
//    friend class boost::serialization::access;

//    template<class Archive>
//    void serialize(Archive & ar, const unsigned int )
//    {
//        ar & m_name;
//        ar & m_abbreviation;
//        ar & m_index;
//        ar & m_mass;
//    }
};

inline void AtomType::setName(string name) {
    m_name = name;
}

inline void AtomType::setAbbreviation(string abbreviation) {
    m_abbreviation = abbreviation;
}

inline void AtomType::setNumber(int id) {
    m_number = id;
}

inline void AtomType::setMass(double mass) {
    m_mass = mass;
}

inline int AtomType::number() const
{
    return m_number;
}

inline string AtomType::name() const
{
    return m_name;
}

inline string AtomType::abbreviation() const
{
    return m_abbreviation;
}

inline double AtomType::effectiveCharge() const
{
    return m_effectiveCharge;
}

inline double AtomType::electronicPolarizability() const
{
    return m_electronicPolarizability;
}

inline void AtomType::setEffectiveCharge(double arg)
{
    m_effectiveCharge = arg;
}

inline void AtomType::setElectronicPolarizability(double arg)
{
    m_electronicPolarizability = arg;
}

inline double AtomType::mass() const
{
    return m_mass;
}

inline int AtomType::index() const
{
    return m_index;
}

inline bool AtomType::operator ==(const AtomType &other) const
{
    return (m_index == other.index());
}

#endif // ATOMICTYPE_H
