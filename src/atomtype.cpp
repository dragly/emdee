#include <src/atomtype.h>

#include <iostream>

AtomType::AtomType() :
    m_name(""),
    m_abbreviation(""),
    m_id(-1),
    m_mass(-1)
{
}

AtomType::AtomType(AtomType::AtomTypeEnum atomTypeEnum)
{
    switch(atomTypeEnum) {
    case HYDROGEN:
        m_name = "Hydrogen";
        m_abbreviation = "H";
        m_id = 1;
        m_mass = -1;
        break;
    case HELIUM:
        m_name = "Helium";
        m_abbreviation = "He";
        m_id = 2;
        m_mass = -1;
        break;
    case OXYGEN:
        m_name = "Oxygen";
        m_abbreviation = "O";
        m_id = 16;
        m_mass = -1;
        break;
    case ARGON:
        m_name = "Argon";
        m_abbreviation = "Ar";
        m_id = 18;
        m_mass = 1; // TODO set to real mass and fix with unit mass
        break;
    default:
        m_name = "Unknown";
        m_abbreviation = "NNN";
        m_id = -1;
        m_mass = -1;
    }
}

AtomType AtomType::argon() {
    return AtomType(ARGON);
}
