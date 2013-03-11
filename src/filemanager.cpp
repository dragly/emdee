#include "filemanager.h"

#include <src/moleculesystem.h>
#include <src/moleculesystemcell.h>
#include <src/atom.h>
#include <src/atomtype.h>
#include <src/integrator/integrator.h>
#include <src/processor.h>

#include <unistd.h>
#include <sys/types.h>
#include <pwd.h>

// System includes
#include <iomanip>
#include <armadillo>
#include <H5Cpp.h>
#include <H5File.h>
//#include <hdf5.h>
#include <boost/filesystem.hpp>

using namespace std;
using namespace arma;
//using namespace H5;

typedef struct s1_t {
    char atomType[3];
    double positionX;
    double positionY;
    double positionZ;
    double velocityX;
    double velocityY;
    double velocityZ;
    double forceX;
    double forceY;
    double forceZ;
    int cellID;
} s1_t;

FileManager::FileManager(MoleculeSystem* system) :
    m_outFileName("/tmp/data*.bin"),
    m_moleculeSystem(system),
    m_unitLength(1),
    m_unitTime(1),
    m_unitEnergy(1),
    m_unitMass(1),
    m_unitTemperature(1)
{
}

bool FileManager::load(string fileName) {
    fileName = parseFileName(fileName);
    if(fileName.find(".xyz") != string::npos) {
        //        return loadXyz(fileName);
        cerr << "HDF5 saving is not yet implemented" << endl;
    } else if(fileName.find(".bin") != string::npos) {
        return loadBinary(fileName);
    } else if(fileName.find(".h5") != string::npos) {
        //        return saveHDF5(step);
        cerr << "HDF5 saving is not yet implemented" << endl;
    }
    return false;
}

bool FileManager::loadBinary(string fileName) {
    cout << "Loading binary file " << fileName << endl;
    ifstream inFile;
    inFile.open(fileName, ios::binary | ios::in);
    if(!inFile.good()) {
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

    //    for(int i = 0; i < 2; i++) {
    //        for(int iDim = 0; iDim < 3; iDim++) {
    //            systemBoundaries[i*3 + iDim] = m_moleculeSystem->boundaries()(i, iDim);
    //        }
    //    }

    inFile.read((char*)&rank, sizeof(int));
    inFile.read((char*)&nProcessors, sizeof(int));

    if(nProcessors > 1) {
        cerr << "Can currently not load more than one processor!" << endl;
        throw new exception();
    }

    inFile.read((char*)&step, sizeof(int));
    inFile.read((char*)&time, sizeof(double));
    inFile.read((char*)&timeStep, sizeof(double));
    inFile.read((char*)&nAtoms, sizeof(int));
    inFile.read((char*)&temperature, sizeof(double));
    inFile.read((char*)&pressure, sizeof(double));
    inFile.read((char*)&averageDisplacement, sizeof(double));
    inFile.read((char*)&averageSquareDisplacement, sizeof(double));
    inFile.read((char*)&systemBoundaries[0], sizeof(double) * 6);

    // Convert back units
    time /= m_unitTime;
    timeStep /= m_unitTime;
    temperature /= m_unitTemperature;
    averageDisplacement /= m_unitLength;
    averageSquareDisplacement /= (m_unitLength * m_unitLength);

    for(int i = 0; i < 2; i++) {
        for(int iDim = 0; iDim < 3; iDim++) {
            systemBoundaries[i*3 + iDim] /= m_unitLength;
        }
    }

    cout << "step: " << step << " time: " << time << " timeStep: " << timeStep << " nAtoms: " << nAtoms << " temperature: " << temperature << endl << "boundaries:" << endl;
    cout << systemBoundaries[0] << " "  << systemBoundaries[1] << " "  << systemBoundaries[2] << " " << endl;
    cout << systemBoundaries[3] << " "  << systemBoundaries[4] << " "  << systemBoundaries[5] << " " << endl;

    // Append step by one
    m_moleculeSystem->setStep(step + 1);
    m_moleculeSystem->integrator()->setTimeStep(timeStep);
    m_moleculeSystem->setTime(time);
    m_moleculeSystem->setBoundaries(systemBoundaries[0], systemBoundaries[3], systemBoundaries[1], systemBoundaries[4], systemBoundaries[2], systemBoundaries[5]);


    vector<Atom*> atoms;
    // Read atom data
    m_moleculeSystem->deleteAtoms();
    cout << "Reading data for " << nAtoms << " atoms" << endl;
//    if(nAtoms > 10000) {
//        cout << "Too many atoms, "
//        exit(9999);
//    }
    for(int i = 0; i < nAtoms; i++) {
        //        cout << "Reading atom " << i << endl;
        Atom* atom = new Atom(AtomType::argon());
        Vector3 position;
        Vector3 velocity;
        Vector3 force;
        Vector3 displacement;
        double potential = 0;
        bool isPositionFixed;
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
        inFile.read((char*)&isPositionFixed, sizeof(bool));
        inFile.read((char*)&id, sizeof(int));

        //        cout << "Atom type is " << atomType << endl;

        position /= m_unitLength;
        displacement /= m_unitLength;
        velocity /= (m_unitLength / m_unitTime);
        force /= (m_unitLength * m_unitMass / (m_unitTime * m_unitTime));

        atom->setPosition(position);
        atom->setVelocity(velocity);
        atom->addForce(force);
        atom->clearDisplacement();
        atom->addDisplacement(displacement);
        atom->addPotential(potential);
        atom->setPositionFixed(isPositionFixed);
        atoms.push_back(atom);
    }

    m_moleculeSystem->addAtoms(atoms);
    inFile.close();
    return true;
}
bool FileManager::saveBinary(int step) {
    stringstream outStepName;
    outStepName << setw(6) << setfill('0') << step;
    string outFileNameLocal = m_outFileName;
    unsigned found = m_outFileName.find_last_of("/\\");
    string outPath = m_outFileName.substr(0,found);
    boost::filesystem::create_directories(outPath);
    size_t starPos = m_outFileName.find("*");
    stringstream processorName;
    processorName << "." << setw(4) << setfill('0') << m_moleculeSystem->processor()->rank();
    if(m_moleculeSystem->processor()->rank() != 0) {
        outFileNameLocal.append(processorName.str());
    }
    ofstream outFile;
    if(starPos != string::npos) {
        outFileNameLocal.replace(starPos, 1, outStepName.str());
        outFile.open(outFileNameLocal, ios::out | ios::binary);
    } else {
        outFile.open(outFileNameLocal, ios::app | ios::out | ios::binary);
    }

//    cout << "Writing to file " << outFileNameLocal << endl;

    // Write header data
    int rank = m_moleculeSystem->processor()->rank();
    int nProcessors = m_moleculeSystem->processor()->nProcessors();
    double time = m_moleculeSystem->time() * m_unitTime;
    double timeStep = m_moleculeSystem->integrator()->timeStep() * m_unitTime;
    int nAtoms = m_moleculeSystem->processor()->nAtoms();
    double temperature = m_moleculeSystem->temperature() * m_unitTemperature;
    double pressure = m_moleculeSystem->pressure() * m_unitMass / (m_unitLength * m_unitTime * m_unitTime);
    double averageDisplacement = m_moleculeSystem->averageDisplacement() * m_unitLength;
    double averageSquareDisplacement = m_moleculeSystem->averageSquareDisplacement() * m_unitLength * m_unitLength;

    double systemBoundaries[6];

//    cout << "Writing temperature of " << temperature << endl;

    for(int i = 0; i < 2; i++) {
        for(int iDim = 0; iDim < 3; iDim++) {
            systemBoundaries[i*3 + iDim] = m_moleculeSystem->boundaries()(i, iDim) * m_unitLength;
        }
    }

    outFile.write((char*)&rank, sizeof(int));
    outFile.write((char*)&nProcessors, sizeof(int));

    outFile.write((char*)&step, sizeof(int));
    outFile.write((char*)&time, sizeof(double));
    outFile.write((char*)&timeStep, sizeof(double));
    outFile.write((char*)&nAtoms, sizeof(int));
    outFile.write((char*)&temperature, sizeof(double));
    outFile.write((char*)&pressure, sizeof(double));
    outFile.write((char*)&averageDisplacement, sizeof(double));
    outFile.write((char*)&averageSquareDisplacement, sizeof(double));
    outFile.write((char*)&systemBoundaries[0], sizeof(double) * 6);

    // Write system boundaries

    // Write atom data
    for(MoleculeSystemCell* cell : m_moleculeSystem->processor()->cells()) {
        for(Atom* atom : cell->atoms()) {
            Vector3 position = atom->position() * m_unitLength;
            Vector3 velocity = atom->velocity() * (m_unitLength / m_unitTime);
            Vector3 force = atom->force() * (m_unitMass * m_unitLength / (m_unitTime * m_unitTime));
            Vector3 displacement = atom->displacement() * m_unitLength;
            double potential = atom->potential() * (m_unitMass * m_unitLength * m_unitLength / (m_unitTime * m_unitTime));
            bool isPositionFixed = atom->isPositionFixed();
            int id = atom->cellID();
            char atomType[3];
            sprintf(atomType, "Ar");
            outFile.write(atomType, sizeof(atomType));
            outFile.write((char*)&position(0), sizeof(double));
            outFile.write((char*)&position(1), sizeof(double));
            outFile.write((char*)&position(2), sizeof(double));
            outFile.write((char*)&velocity(0), sizeof(double));
            outFile.write((char*)&velocity(1), sizeof(double));
            outFile.write((char*)&velocity(2), sizeof(double));
            outFile.write((char*)&force(0), sizeof(double));
            outFile.write((char*)&force(1), sizeof(double));
            outFile.write((char*)&force(2), sizeof(double));
            outFile.write((char*)&displacement(0), sizeof(double));
            outFile.write((char*)&displacement(1), sizeof(double));
            outFile.write((char*)&displacement(2), sizeof(double));
            outFile.write((char*)&potential, sizeof(double));
            outFile.write((char*)&isPositionFixed, sizeof(bool));
            outFile.write((char*)&id, sizeof(int));
        }
    }
    outFile.close();
    return true;
}

bool FileManager::save(int step) {
    // Refreshing to make sure we don't save any wrong cell IDs
    if(m_outFileName.find(".xyz") != string::npos) {
        return saveXyz(step);
    } else if(m_outFileName.find(".bin") != string::npos) {
        return saveBinary(step);
    } else if(m_outFileName.find(".h5") != string::npos) {
        return saveHDF5(step);
    }
    return false;
}

bool FileManager::saveXyz(int step) {
    stringstream outStepName;
    outStepName << setw(6) << setfill('0') << step;
    string outFileNameLocal = m_outFileName;
    size_t starPos = m_outFileName.find("*");
    ofstream outFile;
    if(starPos != string::npos) {
        outFileNameLocal.replace(starPos, 1, outStepName.str());
        outFile.open(outFileNameLocal);
    } else {
        outFile.open(outFileNameLocal, ios_base::app);
    }

    outFile << "Some nice comment" << endl;

    char line[1000];
    for(MoleculeSystemCell* cell : m_moleculeSystem->processor()->cells()) {
        for(Atom* atom : cell->atoms()) {
            Vector3 position = atom->position() * m_unitLength;
            Vector3 velocity = atom->velocity() * (m_unitLength / m_unitTime);
            Vector3 force = atom->force() * (m_unitMass * m_unitLength / (m_unitTime * m_unitTime));
            sprintf(line, " %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %d\n",
                    position(0), position(1), position(2),
                    velocity(0), velocity(1), velocity(2),
                    force(0), force(1), force(2),
                    atom->cellID());
            outFile << atom->type().abbreviation
                    << line;
        }
    }
    outFile.close();
    return true;
}


bool FileManager::saveHDF5(int step) {
    cerr << "saveHDF5 is not yet fully implemented. Missing some parameters and rescaling with units." << endl;
    stringstream outStepName;
    outStepName << setw(6) << setfill('0') << step;
    string outFileNameLocal = m_outFileName;
    size_t starPos = m_outFileName.find("*");
    if(starPos != string::npos) {
        outFileNameLocal.replace(starPos, 1, outStepName.str());
    }
    /*
      * Initialize the data
      */
    // TODO Rescale with units!
    int  i = 0;
    s1_t* s1 = new s1_t[m_moleculeSystem->processor()->nAtoms()];
    for(MoleculeSystemCell* cell : m_moleculeSystem->processor()->cells()) {
        for(Atom* atom : cell->atoms()) {
            sprintf(s1[i].atomType , "Ar");
            s1[i].positionX = atom->position()(0);
            s1[i].positionY = atom->position()(1);
            s1[i].positionZ = atom->position()(2);
            s1[i].velocityX = atom->velocity()(0);
            s1[i].velocityY = atom->velocity()(1);
            s1[i].velocityZ = atom->velocity()(2);
            s1[i].forceX = atom->force()(0);
            s1[i].forceY = atom->force()(1);
            s1[i].forceZ = atom->force()(2);
            s1[i].cellID = atom->cellID();
            i++;
        }
    }

    /*
      * Turn off the auto-printing when failure occurs so that we can
      * handle the errors appropriately
      */
    //        Exception::dontPrint();

    /*
      * Create the data space.
      */
    hsize_t dim[] = {m_moleculeSystem->processor()->nAtoms()};   /* Dataspace dimensions */
    H5::DataSpace space( 1, dim );

    /*
      * Create the file.
      */
    H5::H5File* file = new H5::H5File( outFileNameLocal, H5F_ACC_TRUNC );

    /*
      * Create the memory datatype.
      */
    H5::StrType string_type(H5::PredType::C_S1, 3);
    H5::CompType mtype1( sizeof(s1_t) );
    mtype1.insertMember( "atomType", HOFFSET(s1_t, atomType), string_type);
    mtype1.insertMember( "positionX", HOFFSET(s1_t, positionX), H5::PredType::NATIVE_DOUBLE);
    mtype1.insertMember( "positionY", HOFFSET(s1_t, positionY), H5::PredType::NATIVE_DOUBLE);
    mtype1.insertMember( "positionZ", HOFFSET(s1_t, positionZ), H5::PredType::NATIVE_DOUBLE);
    mtype1.insertMember( "velocityX", HOFFSET(s1_t, velocityX), H5::PredType::NATIVE_DOUBLE);
    mtype1.insertMember( "velocityY", HOFFSET(s1_t, velocityY), H5::PredType::NATIVE_DOUBLE);
    mtype1.insertMember( "velocityZ", HOFFSET(s1_t, velocityZ), H5::PredType::NATIVE_DOUBLE);
    mtype1.insertMember( "forceX", HOFFSET(s1_t, forceX), H5::PredType::NATIVE_DOUBLE);
    mtype1.insertMember( "forceY", HOFFSET(s1_t, forceY), H5::PredType::NATIVE_DOUBLE);
    mtype1.insertMember( "forceZ", HOFFSET(s1_t, forceZ), H5::PredType::NATIVE_DOUBLE);
    mtype1.insertMember( "cellID", HOFFSET(s1_t, cellID), H5::PredType::NATIVE_INT);

    /*
      * Create the dataset.
      */
    H5::DataSet* dataset;
    dataset = new H5::DataSet(file->createDataSet("ArrayOfStructures", mtype1, space));

    /*
      * Write data to the dataset;
      */
    dataset->write( s1, mtype1 );

    /*
      * Release resources
      */
    delete dataset;
    delete file;
    return true;
}

string FileManager::parseFileName(string fileName) {
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

void FileManager::setOutFileName(string fileName)
{
    fileName = parseFileName(fileName);
    // Continue
    m_outFileName = fileName;
}


void FileManager::setUnitLength(double unitLength)
{
    m_unitLength = unitLength;
}

void FileManager::setUnitTime(double unitTime)
{
    m_unitTime = unitTime;
}

void FileManager::setUnitMass(double unitMass)
{
    m_unitMass = unitMass;
}

void FileManager::setUnitTemperature(double unitTemperature) {
    m_unitTemperature = unitTemperature;
}
