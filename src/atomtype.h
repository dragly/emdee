#ifndef ATOMICTYPE_H
#define ATOMICTYPE_H

#include <string>

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
};

#endif // ATOMICTYPE_H
