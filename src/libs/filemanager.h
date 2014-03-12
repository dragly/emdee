#ifndef FILEMANAGER_H
#define FILEMANAGER_H

// Forward
class MoleculeSystem;

// Local

// System
#include <string>
#include <armadillo>
using namespace std;
using namespace arma;

class FileManager
{
public:
    FileManager(MoleculeSystem *system);

    enum FileFormat {
        XyzFormat,
        HDF5Format
    };

    bool load(string fileName);
    bool loadBinary(string fileName);
    bool save(int step);
    bool saveXyz(int step);
//    bool saveHDF5(int step);
    bool saveBinary(int step);

    void setOutFileName(string fileName);

    // Units
    void setUnitLength(double unitLength);
    void setUnitTime(double unitTime);
    void setUnitMass(double unitMass);

    void setUnitTemperature(double unitTemperature);
    void setConfigurationName(string configurationName);
    string parseFileName(string fileName);
    bool setLatestSymlink(int step);
    string headerFileNameFromStep(int step);
    string processorName();
    string lammpsFileNameFromStep(int step);
protected:
    string m_outFileName;

    FileFormat outFileFormat;

    MoleculeSystem *m_moleculeSystem;


    // Units
    double m_unitLength;
    double m_unitTime;
    double m_unitEnergy;
    double m_unitMass;
    double m_unitTemperature;

    string m_configurationName;
private:
    int collectNAtoms();
};

inline void FileManager::setConfigurationName(string configurationName) {
    m_configurationName = configurationName;
}

#endif // FILEMANAGER_H
