#include <src/atomtype.h>

#include <iostream>

AtomType::AtomType() :
    name(""),
    abbreviation(""),
    number(-1)
{
}

AtomType::AtomType(AtomType::AtomTypeEnum atomTypeEnum)
{
    switch(atomTypeEnum) {
    case HYDROGEN:
        name = "Hydrogen";
        abbreviation = "H";
        number = 1;
        mass = -1;
        break;
    case HELIUM:
        name = "Helium";
        abbreviation = "He";
        number = 2;
        mass = -1;
        break;
    case OXYGEN:
        name = "Oxygen";
        abbreviation = "O";
        number = 16;
        mass = -1;
        break;
    case ARGON:
        name = "Argon";
        abbreviation = "Ar";
        number = 18;
        mass = 1;
        break;
    default:
        name = "Unknown";
        abbreviation = "NNN";
        number = -1;
        mass = -1;
    }
}

AtomType AtomType::argon() {
    return AtomType(ARGON);
}
