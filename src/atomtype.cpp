#include <src/atomtype.h>

#include <iostream>

AtomType::AtomType() :
    m_name(""),
    m_abbreviation(""),
    m_number(-1),
    m_mass(-1),
    m_index(-1)
{
}

AtomType::AtomType(int index) :
    m_name(""),
    m_abbreviation(""),
    m_number(-1),
    m_mass(-1),
    m_index(index)
{
}

AtomType::AtomType(AtomType::AtomTypeEnum atomTypeEnum, int index) :
    m_index(index)
{
    switch(atomTypeEnum) {
    case HYDROGEN:
        m_name = "Hydrogen";
        m_abbreviation = "H";
        m_number = 1;
        m_mass = -1;
        break;
    case HELIUM:
        m_name = "Helium";
        m_abbreviation = "He";
        m_number = 2;
        m_mass = -1;
        break;
    case OXYGEN:
        m_name = "Oxygen";
        m_abbreviation = "O";
        m_number = 16;
        m_mass = -1;
        break;
    case ARGON:
        m_name = "Argon";
        m_abbreviation = "Ar";
        m_number = 18;
        m_mass = 1; // TODO set to real mass and fix with unit mass
        break;
    default:
        m_name = "Unknown";
        m_abbreviation = "NNN";
        m_number = -1;
        m_mass = -1;
    }
}

AtomType AtomType::argon() {
    return AtomType(ARGON, 18);
}
