#include <iostream>
#include <fstream>
#include <iomanip>
#include <armadillo>

#include <unistd.h>
#include <sys/types.h>
#include <pwd.h>

using namespace std;
using namespace arma;

string parseFileName(string fileName) {
    // Replace tilde ~ with home directory
    struct passwd *pw = getpwuid(getuid());
    const char *homedir = pw->pw_dir;

    stringstream homeDirName;
    homeDirName << homedir;
    size_t tildePos = fileName.find("~");
    if(tildePos != string::npos) {
        fileName.replace(tildePos, 1, homeDirName.str());
    }
    return fileName;
}

bool convertFile(string fileName) {

    string fileEnding = ".xyz";
    string outFileNameLocal = fileName;
    size_t starPos = fileName.find(".bin");
    ofstream outFile;
    outFileNameLocal.replace(starPos, 4, fileEnding);
    outFile.open(outFileNameLocal);
    //    cout << "Saving xyz data to " << outFileNameLocal  << endl;

    cout << "Loading binary file " << fileName << " and writing to " << outFileNameLocal << endl;
    fileName = parseFileName(fileName);
    ifstream inFile;
    inFile.open(fileName, ios::binary | ios::in);
    if(!inFile.good()) {
        return false;
    }
    // Read header data
    int step;
    double time; // = step * m_moleculeSystem->integrator()->timeStep() * m_unitTime;
    double timeStep; // = m_moleculeSystem->integrator()->timeStep() * m_unitTime;
    int nAtoms; // = m_moleculeSystem->atoms().size();
    double temperature; // = m_moleculeSystem->temperature() * m_unitTemperature;
    double pressure; // = m_moleculeSystem->temperature() * m_unitTemperature;
    double averageDisplacement; // = m_moleculeSystem->averageDisplacement() * m_unitLength;
    double averageSquareDisplacement; // = m_moleculeSystem->averageSquareDisplacement() * m_unitLength * m_unitLength;

    double systemBoundaries[6];

    //    for(int i = 0; i < 2; i++) {
    //        for(int iDim = 0; iDim < 3; iDim++) {
    //            systemBoundaries[i*3 + iDim] = m_moleculeSystem->boundaries()(i, iDim);
    //        }
    //    }

    inFile.read((char*)&step, sizeof(int));
    inFile.read((char*)&time, sizeof(double));
    inFile.read((char*)&timeStep, sizeof(double));
    inFile.read((char*)&nAtoms, sizeof(int));
    inFile.read((char*)&temperature, sizeof(double));
    inFile.read((char*)&pressure, sizeof(double));
    inFile.read((char*)&averageDisplacement, sizeof(double));
    inFile.read((char*)&averageSquareDisplacement, sizeof(double));
    inFile.read((char*)&systemBoundaries[0], sizeof(double) * 6);



    outFile << nAtoms << endl;
    outFile << "Some nice comment" << endl;

//    cout << "Reading data for " << nAtoms << " atoms" << endl;
    if(nAtoms > 100000) {
        exit(9999);
    }
    char line[1000];

    double unitLength = 1e-10;

    for(int i = 0; i < nAtoms; i++) {
//        cout << "Reading atom " << i << endl;
        rowvec position = zeros<rowvec>(3);
        rowvec velocity = zeros<rowvec>(3);
        rowvec force = zeros<rowvec>(3);
        rowvec displacement = zeros<rowvec>(3);
        double potential = 0;
        int id = 0;
        char atomType[3];
        inFile.read(atomType, sizeof(atomType));
        inFile.read((char*)&position(0), sizeof(double));
        inFile.read((char*)&position(1), sizeof(double));
        inFile.read((char*)&position(2), sizeof(double));
        inFile.read((char*)&velocity(0), sizeof(double));
        inFile.read((char*)&velocity(1), sizeof(double));
        inFile.read((char*)&velocity(2), sizeof(double));
        inFile.read((char*)&force(0), sizeof(double));
        inFile.read((char*)&force(1), sizeof(double));
        inFile.read((char*)&force(2), sizeof(double));
        inFile.read((char*)&displacement(0), sizeof(double));
        inFile.read((char*)&displacement(1), sizeof(double));
        inFile.read((char*)&displacement(2), sizeof(double));
        inFile.read((char*)&potential, sizeof(double));
        inFile.read((char*)&id, sizeof(int));

        position /= unitLength;

        sprintf(line, " %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %d\n",
                position(0), position(1), position(2),
                velocity(0), velocity(1), velocity(2),
                force(0), force(1), force(2),
                potential,
                id);
        //                toPrint << atom->type().abbreviation
        //                        << " " << position(0) << " " << position(1) << " " << position(2)
        //                        << " " << velocity(0) << " " << velocity(1) << " " << velocity(2)
        //                        << " " << force(0) << " " << force(1) << " " << force(2)
        //                        << " " << cell->id()
        //                        << "\n";
        outFile << atomType
                << line;
    }

    outFile.close();
    inFile.close();
    return true;
}

int main(int argc, char* argv[])
{
    for(int i = 1; i < argc; i++) {
        convertFile(argv[i]);
    }
    return 0;
}

