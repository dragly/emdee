#ifndef CONFIGURATIONPARSER_H
#define CONFIGURATIONPARSER_H

class MoleculeSystem;

#include <string>

using namespace std;

class ConfigurationParser
{
public:
    ConfigurationParser(MoleculeSystem* moleculeSystem);

    void runConfiguration(string configurationPath);
protected:
    MoleculeSystem* m_moleculeSystem;
};

#endif // CONFIGURATIONPARSER_H
