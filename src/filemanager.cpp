#include "filemanager.h"

#include <src/moleculesystem.h>
#include <src/moleculesystemcell.h>
#include <src/atom.h>
#include <src/molecule.h>
#include <src/atomtype.h>
#include <src/integrator/integrator.h>

#include <unistd.h>
#include <sys/types.h>
#include <pwd.h>

// System includes
#include <iomanip>
#include <armadillo>
#include <H5Cpp.h>
#include <H5File.h>
//#include <hdf5.h>

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
    } else if(fileName.find(".bin") != string::npos) {
        return loadBinary(fileName);
    } else if(fileName.find(".h5") != string::npos) {
        //        return saveHDF5(step);
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
    int step;
    double time; // = step * m_moleculeSystem->integrator()->timeStep() * m_unitTime;
    double timeStep; // = m_moleculeSystem->integrator()->timeStep() * m_unitTime;
    int nAtoms; // = m_moleculeSystem->atoms().size();
    double temperature; // = m_moleculeSystem->temperature() * m_unitTemperature;
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

    cout << step << " " << time << " " << timeStep << " " << nAtoms << " " << temperature << " " << endl;
    cout << systemBoundaries[0] << " "  << systemBoundaries[1] << " "  << systemBoundaries[2] << " " << endl;
    cout << systemBoundaries[3] << " "  << systemBoundaries[4] << " "  << systemBoundaries[5] << " " << endl;

    // Append step by one
    m_moleculeSystem->setStep(step + 1);
    m_moleculeSystem->integrator()->setTimeStep(timeStep);
    m_moleculeSystem->setBoundaries(systemBoundaries[0], systemBoundaries[3], systemBoundaries[1], systemBoundaries[4], systemBoundaries[2], systemBoundaries[5]);


    vector<Molecule*> molecules;
    // Read atom data
    m_moleculeSystem->deleteMoleculesAndAtoms();
    cout << "Reading data for " << nAtoms << " atoms" << endl;
    if(nAtoms > 10000) {
        exit(9999);
    }
    for(int i = 0; i < nAtoms; i++) {
//        cout << "Reading atom " << i << endl;
        Molecule* molecule = new Molecule();
        Atom* atom = new Atom(molecule, AtomType::argon());
        molecule->addAtom(atom);
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

//        cout << "Atom type is " << atomType << endl;

        position /= m_unitLength;
        displacement /= m_unitLength;
        velocity /= (m_unitLength / m_unitTime);
        force /= (m_unitLength * m_unitMass / (m_unitTime * m_unitTime));

        molecule->setPosition(position);
        molecule->setVelocity(velocity);
        atom->addForce(force);
        molecule->clearDisplacement();
        molecule->addDisplacement(displacement);
        atom->addPotential(potential);
        molecules.push_back(molecule);
    }

    m_moleculeSystem->addMolecules(molecules);
    inFile.close();
    return true;
}

bool FileManager::save(int step) {
    // Refreshing to make sure we don't save any wrong cell IDs
    //    refreshCellContents();

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
    //    cout << "Saving xyz data to " << outFileNameLocal  << endl;

    outFile << m_moleculeSystem->molecules().size() << endl;
    outFile << "Some nice comment" << endl;

    char line[1000];
    for(MoleculeSystemCell* cell : m_moleculeSystem->cells()) {
        for(Molecule* molecule : cell->molecules()) {
            for(Atom* atom : molecule->atoms()) {
                rowvec position = atom->position() * m_unitLength;
                rowvec velocity = atom->velocity() * (m_unitLength / m_unitTime);
                rowvec force = atom->force() * (m_unitMass * m_unitLength / (m_unitTime * m_unitTime));
                sprintf(line, " %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %d\n",
                        position(0), position(1), position(2),
                        velocity(0), velocity(1), velocity(2),
                        force(0), force(1), force(2),
                        cell->id());
                //                toPrint << atom->type().abbreviation
                //                        << " " << position(0) << " " << position(1) << " " << position(2)
                //                        << " " << velocity(0) << " " << velocity(1) << " " << velocity(2)
                //                        << " " << force(0) << " " << force(1) << " " << force(2)
                //                        << " " << cell->id()
                //                        << "\n";
                outFile << atom->type().abbreviation
                        << line;
            }
        }
    }
    outFile.close();
    return true;
}

bool FileManager::saveBinary(int step) {
    stringstream outStepName;
    outStepName << setw(6) << setfill('0') << step;
    string outFileNameLocal = m_outFileName;
    size_t starPos = m_outFileName.find("*");
    ofstream outFile;
    if(starPos != string::npos) {
        outFileNameLocal.replace(starPos, 1, outStepName.str());
        outFile.open(outFileNameLocal, ios::out | ios::binary);
    } else {
        outFile.open(outFileNameLocal, ios::app | ios::out | ios::binary);
    }
    //    cout << "Saving binary data to " << outFileNameLocal  << endl;

    //    outFile << m_molecules.size() << endl;
    //    outFile << "Some nice comment" << endl;

    //    char line[1000];
    // Write header data
    double time = step * m_moleculeSystem->integrator()->timeStep() * m_unitTime;
    double timeStep = m_moleculeSystem->integrator()->timeStep() * m_unitTime;
    int nAtoms = m_moleculeSystem->atoms().size();
    double temperature = m_moleculeSystem->temperature() * m_unitTemperature;
    double averageDisplacement = m_moleculeSystem->averageDisplacement() * m_unitLength;
    double averageSquareDisplacement = m_moleculeSystem->averageSquareDisplacement() * m_unitLength * m_unitLength;

    double systemBoundaries[6];

    for(int i = 0; i < 2; i++) {
        for(int iDim = 0; iDim < 3; iDim++) {
            systemBoundaries[i*3 + iDim] = m_moleculeSystem->boundaries()(i, iDim) * m_unitLength;
        }
    }

    outFile.write((char*)&step, sizeof(int));
    outFile.write((char*)&time, sizeof(double));
    outFile.write((char*)&timeStep, sizeof(double));
    outFile.write((char*)&nAtoms, sizeof(int));
    outFile.write((char*)&temperature, sizeof(double));
    outFile.write((char*)&averageDisplacement, sizeof(double));
    outFile.write((char*)&averageSquareDisplacement, sizeof(double));
    outFile.write((char*)&systemBoundaries[0], sizeof(double) * 6);

    // Write system boundaries

    // Write atom data
    for(MoleculeSystemCell* cell : m_moleculeSystem->cells()) {
        for(Molecule* molecule : cell->molecules()) {
            for(Atom* atom : molecule->atoms()) {
                rowvec position = atom->position() * m_unitLength;
                rowvec velocity = atom->velocity() * (m_unitLength / m_unitTime);
                rowvec force = atom->force() * (m_unitMass * m_unitLength / (m_unitTime * m_unitTime));
                rowvec displacement = molecule->displacement() * m_unitLength;
                double potential = atom->potential() * (m_unitMass * m_unitLength * m_unitLength / (m_unitTime * m_unitTime));
                //                sprintf(line, " %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %d\n",
                //                        position(0), position(1), position(2),
                //                        velo0city(0), velocity(1), velocity(2),
                //                        force(0), force(1), force(2),
                //                        cell->id());
                int id = cell->id();
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
                outFile.write((char*)&id, sizeof(int));
            }
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
    s1_t* s1 = new s1_t[m_moleculeSystem->atoms().size()];
    for(MoleculeSystemCell* cell : m_moleculeSystem->cells()) {
        for(Molecule* molecule : cell->molecules()) {
            for(Atom* atom : molecule->atoms()) {
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
                s1[i].cellID = cell->id();
                i++;
            }
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
    hsize_t dim[] = {m_moleculeSystem->atoms().size()};   /* Dataspace dimensions */
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
