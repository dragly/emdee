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
int readNAtoms(string fileName) {
    cout << "Loading binary file " << fileName << " to count atoms" << endl;
    ifstream inFile;
    inFile.open(fileName, ios::binary | ios::in);
    if(!inFile.good()) {
        cout << "File is no good!" << endl;
        return false;
    }
    // Read header data
    int rank;
    int nProcessors;

    int step;
    double time; // = step * m_moleculeSystem->integrator()->timeStep() * m_unitTime;
    double timeStep; // = m_moleculeSystem->integrator()->timeStep() * m_unitTime;
    int nAtoms; // = m_moleculeSystem->atoms().size();
    double temperature; // = m_moleculeSystem->temperature() * m_unitTemperature;
    double pressure; // = m_moleculeSystem->temperature() * m_unitTemperature;
    double averageDisplacement; // = m_moleculeSystem->averageDisplacement() * m_unitLength;
    double averageSquareDisplacement; // = m_moleculeSystem->averageSquareDisplacement() * m_unitLength * m_unitLength;

    double systemBoundaries[6];

    inFile.read((char*)&rank, sizeof(int));
    inFile.read((char*)&nProcessors, sizeof(int));

//    cout << "Reading file for rank " << rank << " of " << nProcessors << endl;

    inFile.read((char*)&step, sizeof(int));
    inFile.read((char*)&time, sizeof(double));
    inFile.read((char*)&timeStep, sizeof(double));
    inFile.read((char*)&nAtoms, sizeof(int));
    inFile.read((char*)&temperature, sizeof(double));
    inFile.read((char*)&pressure, sizeof(double));
    inFile.read((char*)&averageDisplacement, sizeof(double));
    inFile.read((char*)&averageSquareDisplacement, sizeof(double));
    inFile.read((char*)&systemBoundaries[0], sizeof(double) * 6);

//    cout << "Found " << nAtoms << " atoms." << endl;
//    cout << "Found " << temperature << " temperature." << endl;

    inFile.close();
    return nAtoms;
}

bool readWriteFile(string fileName, ofstream &outFile, bool isFirst) {
    cout << "Loading binary file " << fileName << endl;
    ifstream headerFile;
    headerFile.open(fileName, ios::binary | ios::in);
    if(!headerFile.good()) {
        cout << "Header file " << fileName << " no good...";
        return false;
    }
    // Read header data
    int rank;
    int nProcessors;
    int step;
    double time; // = step * m_moleculeSystem->integrator()->timeStep() * m_unitTime;
    double timeStep; // = m_moleculeSystem->integrator()->timeStep() * m_unitTime;
    int nAtoms; // = m_moleculeSystem->atoms().size();
    double temperature; // = m_moleculeSystem->temperature() * m_unitTemperature;
    double pressure;
    double averageDisplacement; // = m_moleculeSystem->averageDisplacement() * m_unitLength;
    double averageSquareDisplacement; // = m_moleculeSystem->averageSquareDisplacement() * m_unitLength * m_unitLength;

    double systemBoundaries[6];
    double kineticEnergyTotal;
    double potentialEnergyTotal;

    //    for(int i = 0; i < 2; i++) {
    //        for(int iDim = 0; iDim < 3; iDim++) {
    //            systemBoundaries[i*3 + iDim] = m_moleculeSystem->boundaries()(i, iDim);
    //        }
    //    }

    headerFile.read((char*)&rank, sizeof(int));
    headerFile.read((char*)&nProcessors, sizeof(int));
    headerFile.close();

    for(int processor = 0; processor < nProcessors; processor++) {
        ifstream lammpsFile;
        string subFileName = fileName;
        string lammpsFileName = fileName;
        lammpsFileName.replace(lammpsFileName.find(".bin"), 4, ".lmp");
        stringstream processorName;
        processorName << "." << setw(4) << setfill('0') << processor;
        if(processor > 0) {
            lammpsFileName.append(processorName.str());
            subFileName.append(processorName.str());
        }
        headerFile.open(subFileName, ios::binary | ios::in);
        lammpsFile.open(lammpsFileName, ios::binary | ios::in);
        if(!lammpsFile.good() || !headerFile.good()) {
            cout << "File " << subFileName << " or " << lammpsFileName << " no good...";
            return false;
        }

        if(processor == 0) {
            int nAtomsTotal;
            double shear;
            int nChunks;
            int chunkLength;
            int nColumns;
            lammpsFile.read((char*)&step, sizeof(int));
            lammpsFile.read((char*)&nAtomsTotal, sizeof(int));
            lammpsFile.read((char*)&systemBoundaries[0], sizeof(double));
            lammpsFile.read((char*)&systemBoundaries[3], sizeof(double));
            lammpsFile.read((char*)&systemBoundaries[1], sizeof(double));
            lammpsFile.read((char*)&systemBoundaries[4], sizeof(double));
            lammpsFile.read((char*)&systemBoundaries[2], sizeof(double));
            lammpsFile.read((char*)&systemBoundaries[5], sizeof(double));
            lammpsFile.read((char*)&shear, sizeof(double));
            lammpsFile.read((char*)&shear, sizeof(double));
            lammpsFile.read((char*)&shear, sizeof(double));
            lammpsFile.read((char*)&nColumns, sizeof(int));
            lammpsFile.read((char*)&nChunks, sizeof(int));
            lammpsFile.read((char*)&chunkLength, sizeof(int));
            outFile << nAtomsTotal << endl;
            outFile << "Comment" << endl;
        }

        headerFile.read((char*)&rank, sizeof(int));
        headerFile.read((char*)&nProcessors, sizeof(int));
        headerFile.read((char*)&step, sizeof(int));
        headerFile.read((char*)&time, sizeof(double));
        headerFile.read((char*)&timeStep, sizeof(double));
        headerFile.read((char*)&nAtoms, sizeof(int));
        headerFile.read((char*)&temperature, sizeof(double));
        headerFile.read((char*)&pressure, sizeof(double));
        headerFile.read((char*)&averageDisplacement, sizeof(double));
        headerFile.read((char*)&averageSquareDisplacement, sizeof(double));
        headerFile.read((char*)&systemBoundaries[0], sizeof(double) * 6);
        headerFile.read((char*)&kineticEnergyTotal, sizeof(double));
        headerFile.read((char*)&potentialEnergyTotal, sizeof(double));

        // Convert back units
//        time /= m_unitTime;
//        timeStep /= m_unitTime;
//        temperature /= m_unitTemperature;
//        averageDisplacement /= m_unitLength;
//        averageSquareDisplacement /= (m_unitLength * m_unitLength);

//        for(int i = 0; i < 2; i++) {
//            for(int iDim = 0; iDim < 3; iDim++) {
//                systemBoundaries[i*3 + iDim] /= m_unitLength;
//            }
//        }

        cout << "step: " << step << " time: " << time << " timeStep: " << timeStep << " nAtoms: " << nAtoms << " temperature: " << temperature << endl << "boundaries:" << endl;
        cout << systemBoundaries[0] << " "  << systemBoundaries[1] << " "  << systemBoundaries[2] << " " << endl;
        cout << systemBoundaries[3] << " "  << systemBoundaries[4] << " "  << systemBoundaries[5] << " " << endl;

        // Append step by one
//        m_moleculeSystem->setStep(step + 1);
//        m_moleculeSystem->integrator()->setTimeStep(timeStep);
//        m_moleculeSystem->setTime(time);
//        m_moleculeSystem->setBoundaries(systemBoundaries[0], systemBoundaries[3], systemBoundaries[1], systemBoundaries[4], systemBoundaries[2], systemBoundaries[5]);


//        vector<Atom*> atoms;
        // Read atom data
        //    m_moleculeSystem->deleteAtoms();
        cout << "Reading data for " << nAtoms << " atoms" << endl;
        //    if(nAtoms > 10000) {
        //        cout << "Too many atoms, "
        //        exit(9999);
        //    }
        for(int i = 0; i < nAtoms; i++) {
            //        cout << "Reading atom " << i << endl;
//            Atom* atom = new Atom(AtomType::argon());
            rowvec position = zeros<rowvec>(3);
            rowvec velocity = zeros<rowvec>(3);
            rowvec force = zeros<rowvec>(3);
            rowvec displacement = zeros<rowvec>(3);;
            double potential = 0;
            double isPositionFixed;
            double id = 0;
            double atomType = 0;

            lammpsFile.read((char*)&atomType, sizeof(double)); // Why is the atomType a double? Because this is the standard in LAMMPS files, which are now used as output
            lammpsFile.read((char*)&position(0), sizeof(double));
            lammpsFile.read((char*)&position(1), sizeof(double));
            lammpsFile.read((char*)&position(2), sizeof(double));
            lammpsFile.read((char*)&velocity(0), sizeof(double));
            lammpsFile.read((char*)&velocity(1), sizeof(double));
            lammpsFile.read((char*)&velocity(2), sizeof(double));
            lammpsFile.read((char*)&force(0), sizeof(double));
            lammpsFile.read((char*)&force(1), sizeof(double));
            lammpsFile.read((char*)&force(2), sizeof(double));
            lammpsFile.read((char*)&displacement(0), sizeof(double));
            lammpsFile.read((char*)&displacement(1), sizeof(double));
            lammpsFile.read((char*)&displacement(2), sizeof(double));
            lammpsFile.read((char*)&potential, sizeof(double));
            lammpsFile.read((char*)&isPositionFixed, sizeof(double));
            lammpsFile.read((char*)&id, sizeof(double));

            cout << isPositionFixed << endl;

            if((int)isPositionFixed > 0) {
                outFile << "Ar ";
            } else {
                outFile << "Ne ";
            }
            outFile << position(0) << " ";
            outFile << position(1) << " ";
            outFile << position(2) << " ";
            outFile << velocity(0) << " ";
            outFile << velocity(1) << " ";
            outFile << velocity(2) << " ";
            outFile << force(0) << " ";
            outFile << force(1) << " ";
            outFile << force(2) << " ";
            outFile << displacement(0) << " ";
            outFile << displacement(1) << " ";
            outFile << displacement(2) << " ";
            outFile << potential << " ";
            outFile << (int)isPositionFixed << " ";
            outFile << (int)id << " ";
            outFile << endl;

            //        cout << "Atom type is " << atomType << endl;

//            position /= (m_unitLength * 1e10);
//            displacement /= m_unitLength;
//            velocity /= (m_unitLength / m_unitTime);
//            force /= (m_unitLength * m_unitMass / (m_unitTime * m_unitTime));

//            atom->setPosition(position);
//            atom->setVelocity(velocity);
//            atom->addForce(force);
//            atom->clearDisplacement();
//            atom->addDisplacement(displacement);
//            atom->addPotential(potential);
//            atom->setPositionFixed(isPositionFixed);
//            atoms.push_back(atom);
        }

//        m_moleculeSystem->addAtoms(atoms);
        headerFile.close();
        lammpsFile.close();
    }
    return true;
}

bool convertFile(string fileName) {

    string fileEnding = ".xyz";
    string outFileNameLocal = fileName;
    size_t starPos = fileName.find(".bin");

    ofstream outFile;
    outFileNameLocal.replace(starPos, 4, fileEnding);
    //    cout << "Saving xyz data to " << outFileNameLocal  << endl;

    fileName = parseFileName(fileName);
    outFileNameLocal = parseFileName(outFileNameLocal);

    outFile.open(outFileNameLocal);

    readWriteFile(fileName, outFile, true);

    outFile.close();
    return true;
}

int main(int argc, char* argv[])
{
    for(int i = 1; i < argc; i++) {
        convertFile(argv[i]);
    }
    return 0;
}

