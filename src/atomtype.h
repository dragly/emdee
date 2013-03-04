#ifndef ATOMICTYPE_H
#define ATOMICTYPE_H

#include <string>

#include <boost/mpi.hpp>
#include <boost/serialization/string.hpp>

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
    AtomType(AtomTypeEnum atomTypeEnum);
    string name;
    string abbreviation;
    int number;
    double mass;
    static AtomType argon();

private:
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive & ar, const unsigned int )
    {
        ar & name;
        ar & abbreviation;
        ar & number;
        ar & mass;
    }
};

#endif // ATOMICTYPE_H
