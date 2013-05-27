#include <iostream>

using namespace std;

#include <src/moleculesystem.h>
#include <src/atom.h>>
#include <src/math/vector3.h>
#include <fstream>
#include <armadillo>
#include <iomanip>
#include <pwd.h>

using namespace arma;
using namespace std;

int main(int argc, char* argv[])
{
    if(argc < 2) {
        cout << "EROR! Too few arguments..." << endl;
        cout << "Usage: bin2xyz inputfile1 [, inputfile2, ...]" << endl;
        exit(0);
    }

    for(int iFile = 0; iFile < argc; iFile++) {
        string fileName = string(argv[iFile]);
        // Replace tilde ~ with home directory
        struct passwd *pw = getpwuid(getuid());
        const char *homedir = pw->pw_dir;

        stringstream homeDirName;
        homeDirName << homedir;
        size_t tildePos = fileName.find("~");
        if(tildePos != string::npos) {
            fileName.replace(tildePos, 1, homeDirName.str());
        }

        MoleculeSystem system;
        vector<AtomType> particleTypes;
        AtomType silicon;
        silicon.setAbbreviation("Si");
        silicon.setName("Silicon");
        silicon.setNumber(14);
        particleTypes.push_back(silicon);
        AtomType oxygen;
        oxygen.setAbbreviation("O");
        oxygen.setName("Oxygen");
        oxygen.setNumber(8);
        particleTypes.push_back(oxygen);
        AtomType argon;
        argon.setAbbreviation("Ar");
        argon.setName("Argon");
        argon.setNumber(18);
        particleTypes.push_back(argon);
        system.setParticleTypes(particleTypes);

        cout << "Loading binary file " << fileName << endl;

        system.load(argv[1]);

        ofstream xyzFile;
        string xyzFileName = fileName;
        size_t suffixPos = xyzFileName.find(".bin");
        if(suffixPos != string::npos) {
            xyzFileName.replace(suffixPos, 4, ".xyz");
        }
        xyzFile.open(xyzFileName, ios::out);

        xyzFile << system.atoms().size() << endl;
        xyzFile << "Output from bin2xyz" << endl;
        for(uint i = 0; i < system.atoms().size(); i++) {
            Atom* atom1 = system.atoms()[i];
            Vector3 position = atom1->position() * 1e10;
            xyzFile << atom1->type().abbreviation() << setprecision(15) << " " << position(0) << " " << position(1) << " " << position(2) << endl;
        }

        xyzFile.close();
    }
    return 0;
}

