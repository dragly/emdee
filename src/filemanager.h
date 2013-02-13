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

    void load(string fileName);
    bool save(int step);
    bool saveXyz(int step);
    bool saveHDF5(int step);
    bool saveBinary(int step);

    void setOutFileName(string fileName);

private:
    string m_outFileName;

    FileFormat outFileFormat;

    MoleculeSystem *m_moleculeSystem;

};

#endif // FILEMANAGER_H
