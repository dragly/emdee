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
#include <unistd.h>
#include <iomanip>
#include <armadillo>
//#include <H5Cpp.h>
//#include <H5File.h>
//#include <hdf5.h>
#ifdef USE_MPI
#include <boost/filesystem.hpp>
#include <boost/mpi.hpp>
namespace mpi = boost::mpi;
#endif

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
    m_unitLength(1e-10),
    m_unitTime(1e-15),
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
        //    m_moleculeSystem->deleteAtoms();
        cout << "Reading data for " << nAtoms << " atoms" << endl;
        //    if(nAtoms > 10000) {
        //        cout << "Too many atoms, "
        //        exit(9999);
        //    }
        const unordered_map<int,AtomType>& particleTypes = m_moleculeSystem->particleTypesById();
        for(int i = 0; i < nAtoms; i++) {
            //        cout << "Reading atom " << i << endl;
            double atomID = -1;
            Vector3 position;
            Vector3 velocity;
            Vector3 force;
            Vector3 displacement;
            double potential = 0;
            double isPositionFixed;
            double cellID = 0;
            double atomTypeId = 0;

             // Why is the atomType a double? Because this is the standard in LAMMPS files, which are now used as output
            lammpsFile.read((char*)&atomID, sizeof(double));
            lammpsFile.read((char*)&atomTypeId, sizeof(double));
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
            lammpsFile.read((char*)&cellID, sizeof(double));

            //        cout << "Atom type is " << atomType << endl;

            position /= (m_unitLength * 1e10);
            displacement /= m_unitLength;
            velocity /= (m_unitLength / m_unitTime);
            force /= (m_unitLength * m_unitMass / (m_unitTime * m_unitTime));

            Atom* atom;
            if(particleTypes.find(atomTypeId) == particleTypes.end()) {
                AtomType atomType;
                atomType.setNumber(atomTypeId);
                atom = new Atom(atomType);
            } else {
                const AtomType& atomType = particleTypes.at(static_cast<int>(atomTypeId));
                atom = new Atom(atomType);
            }

            atom->setID(atomID);
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
        headerFile.close();
        lammpsFile.close();
    }
    return true;
}

bool FileManager::setLatestSymlink(int step) {
    string headerFileName = headerFileNameFromStep(step);
    string lammpsFileName = lammpsFileNameFromStep(step);

    // Create symlink
    string headerSymlinkFileName = headerFileName.substr(0, headerFileName.find_last_of("/") + 1);
    headerSymlinkFileName.append(m_configurationName);
    headerSymlinkFileName.replace(headerSymlinkFileName.find(".cfg"), 4, ".bin");
    string lammpsSymlinkFileName = headerSymlinkFileName;
    lammpsSymlinkFileName.replace(lammpsSymlinkFileName.find(".bin"), 4, ".lmp");

    headerSymlinkFileName.append(processorName());
    lammpsSymlinkFileName.append(processorName());

    unlink(headerSymlinkFileName.c_str());
    unlink(lammpsSymlinkFileName.c_str());

    string relativeHeaderFileName = headerFileName.substr(headerFileName.find_last_of("/") + 1, -1);
    string relativeLammpsFileName = lammpsFileName.substr(lammpsFileName.find_last_of("/") + 1, -1);

    // Symlink returns 0 on success
    if(symlink(relativeHeaderFileName.c_str(), headerSymlinkFileName.c_str()) || symlink(relativeLammpsFileName.c_str(), lammpsSymlinkFileName.c_str())) {
        cerr << "Could not create symlink" << endl;
        return false;
    }
    return true;
}

string FileManager::processorName() {
    stringstream processorNameLocal;
    if( m_moleculeSystem->processor()->rank() != 0) {
        processorNameLocal << "." << setw(4) << setfill('0') << m_moleculeSystem->processor()->rank();
    }
    return processorNameLocal.str();
}

string FileManager::headerFileNameFromStep(int step) {
    stringstream outStepName;
    outStepName << setw(6) << setfill('0') << step;
    string headerFileName = m_outFileName;
    unsigned found = m_outFileName.find_last_of("/\\");
    string outPath = m_outFileName.substr(0,found);
    boost::filesystem::create_directories(outPath);
    size_t starPos = m_outFileName.find("*");
    headerFileName.append(processorName());
    headerFileName.replace(starPos, 1, outStepName.str());
    return headerFileName;
}

string FileManager::lammpsFileNameFromStep(int step) {
    string headerFileName = headerFileNameFromStep(step);
    string lammpsFileName = headerFileName;
    lammpsFileName.replace(lammpsFileName.find(".bin"), 4, ".lmp");
    return lammpsFileName;
}

bool FileManager::saveBinary(int step) {
    mpi::communicator world;
    ofstream headerFile;
    ofstream lammpsFile;
    string headerFileName = headerFileNameFromStep(step);
    string lammpsFileName = lammpsFileNameFromStep(step);
    // Open files for writing
    headerFile.open(headerFileName, ios::out | ios::binary);
    lammpsFile.open(lammpsFileName, ios::out | ios::binary);

    cout << "Writing to file " << headerFileName << endl;

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
    double kineticEnergyTotal = m_moleculeSystem->kineticEnergyTotal() * m_unitLength * m_unitLength * m_unitMass / (m_unitTime * m_unitTime);
    double potentialEnergyTotal = m_moleculeSystem->potentialEnergyTotal() * m_unitLength * m_unitLength * m_unitMass / (m_unitTime * m_unitTime);

    double systemBoundaries[6];

    //    cout << "Writing temperature of " << temperature << endl;

    for(int i = 0; i < 2; i++) {
        for(int iDim = 0; iDim < 3; iDim++) {
            systemBoundaries[i*3 + iDim] = m_moleculeSystem->boundaries()(i, iDim) * m_unitLength;
        }
    }

    // Write to header file
    headerFile.write((char*)&rank, sizeof(int));
    headerFile.write((char*)&nProcessors, sizeof(int));

    headerFile.write((char*)&step, sizeof(int));
    headerFile.write((char*)&time, sizeof(double));
    headerFile.write((char*)&timeStep, sizeof(double));
    headerFile.write((char*)&nAtoms, sizeof(int));
    headerFile.write((char*)&temperature, sizeof(double));
    headerFile.write((char*)&pressure, sizeof(double));
    headerFile.write((char*)&averageDisplacement, sizeof(double));
    headerFile.write((char*)&averageSquareDisplacement, sizeof(double));
    headerFile.write((char*)&systemBoundaries[0], sizeof(double) * 6);
    headerFile.write((char*)&kineticEnergyTotal, sizeof(double));
    headerFile.write((char*)&potentialEnergyTotal, sizeof(double));

    // Write to LAMMPS type file
    int nColumns = 17;
    int nChunks = 1;
    double shear = 0.0;

    if(m_moleculeSystem->processor()->rank() == 0) {
        int nAtomsTotal = nAtoms;
        for(int i = 1; i < m_moleculeSystem->processor()->nProcessors(); i++) {
            int nOtherAtoms = 0;
            world.recv(i, 141, nOtherAtoms);
            nAtomsTotal += nOtherAtoms;
        }

        for(int i = 0; i < 6; i++) {
            systemBoundaries[i] /= 1e-10; // LAMMPS wants lengths in Ångstrøm
        }

//        cout << "The total number of atoms written to file will be " << nAtomsTotal << endl;
        int chunkLength = nAtomsTotal * nColumns;
        lammpsFile.write((char*)&step, sizeof(int));
        lammpsFile.write((char*)&nAtomsTotal, sizeof(int));
        lammpsFile.write((char*)&systemBoundaries[0], sizeof(double));
        lammpsFile.write((char*)&systemBoundaries[3], sizeof(double));
        lammpsFile.write((char*)&systemBoundaries[1], sizeof(double));
        lammpsFile.write((char*)&systemBoundaries[4], sizeof(double));
        lammpsFile.write((char*)&systemBoundaries[2], sizeof(double));
        lammpsFile.write((char*)&systemBoundaries[5], sizeof(double));
        lammpsFile.write((char*)&shear, sizeof(double));
        lammpsFile.write((char*)&shear, sizeof(double));
        lammpsFile.write((char*)&shear, sizeof(double));
        lammpsFile.write((char*)&nColumns, sizeof(int));
        lammpsFile.write((char*)&nChunks, sizeof(int));
        lammpsFile.write((char*)&chunkLength, sizeof(int));
    } else {
        world.send(0, 141, nAtoms);
    }

    // Write atom data
    for(MoleculeSystemCell* cell : m_moleculeSystem->processor()->cells()) {
        for(Atom* atom : cell->atoms()) {
            double atomID = atom->id();
            Vector3 position = atom->position() * m_unitLength;
            Vector3 velocity = atom->velocity() * (m_unitLength / m_unitTime);
            Vector3 force = atom->force() * (m_unitMass * m_unitLength / (m_unitTime * m_unitTime));
            Vector3 displacement = atom->displacement() * m_unitLength;
            double potential = atom->potential() * (m_unitMass * m_unitLength * m_unitLength / (m_unitTime * m_unitTime));
            double isPositionFixed = (double)atom->isPositionFixed();
            double cellID = (double)atom->cellID();
            double atomType = atom->type().number();

            position /= 1e-10; // LAMMPS wants lengths in Ångstrøm

            lammpsFile.write((char*)&atomID, sizeof(double));
            lammpsFile.write((char*)&atomType, sizeof(double));
            lammpsFile.write((char*)&position(0), sizeof(double));
            lammpsFile.write((char*)&position(1), sizeof(double));
            lammpsFile.write((char*)&position(2), sizeof(double));
            lammpsFile.write((char*)&velocity(0), sizeof(double));
            lammpsFile.write((char*)&velocity(1), sizeof(double));
            lammpsFile.write((char*)&velocity(2), sizeof(double));
            lammpsFile.write((char*)&force(0), sizeof(double));
            lammpsFile.write((char*)&force(1), sizeof(double));
            lammpsFile.write((char*)&force(2), sizeof(double));
            lammpsFile.write((char*)&displacement(0), sizeof(double));
            lammpsFile.write((char*)&displacement(1), sizeof(double));
            lammpsFile.write((char*)&displacement(2), sizeof(double));
            lammpsFile.write((char*)&potential, sizeof(double));
            lammpsFile.write((char*)&isPositionFixed, sizeof(double));
            lammpsFile.write((char*)&cellID, sizeof(double));
        }
    }
    headerFile.close();
    lammpsFile.close();
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
            outFile << atom->type().abbreviation()
                    << line;
        }
    }
    outFile.close();
    return true;
}


//bool FileManager::saveHDF5(int step) {
//    cerr << "saveHDF5 is not yet fully implemented. Missing some parameters and rescaling with units." << endl;
//    stringstream outStepName;
//    outStepName << setw(6) << setfill('0') << step;
//    string outFileNameLocal = m_outFileName;
//    size_t starPos = m_outFileName.find("*");
//    if(starPos != string::npos) {
//        outFileNameLocal.replace(starPos, 1, outStepName.str());
//    }
//    /*
//      * Initialize the data
//      */
//    // TODO Rescale with units!
//    int  i = 0;
//    s1_t* s1 = new s1_t[m_moleculeSystem->processor()->nAtoms()];
//    for(MoleculeSystemCell* cell : m_moleculeSystem->processor()->cells()) {
//        for(Atom* atom : cell->atoms()) {
//            sprintf(s1[i].atomType , "Ar");
//            s1[i].positionX = atom->position()(0);
//            s1[i].positionY = atom->position()(1);
//            s1[i].positionZ = atom->position()(2);
//            s1[i].velocityX = atom->velocity()(0);
//            s1[i].velocityY = atom->velocity()(1);
//            s1[i].velocityZ = atom->velocity()(2);
//            s1[i].forceX = atom->force()(0);
//            s1[i].forceY = atom->force()(1);
//            s1[i].forceZ = atom->force()(2);
//            s1[i].cellID = atom->cellID();
//            i++;
//        }
//    }

//    /*
//      * Turn off the auto-printing when failure occurs so that we can
//      * handle the errors appropriately
//      */
//    //        Exception::dontPrint();

//    /*
//      * Create the data space.
//      */
//    hsize_t dim[] = {m_moleculeSystem->processor()->nAtoms()};   /* Dataspace dimensions */
//    H5::DataSpace space( 1, dim );

//    /*
//      * Create the file.
//      */
//    H5::H5File* file = new H5::H5File( outFileNameLocal, H5F_ACC_TRUNC );

//    /*
//      * Create the memory datatype.
//      */
//    H5::StrType string_type(H5::PredType::C_S1, 3);
//    H5::CompType mtype1( sizeof(s1_t) );
//    mtype1.insertMember( "atomType", HOFFSET(s1_t, atomType), string_type);
//    mtype1.insertMember( "positionX", HOFFSET(s1_t, positionX), H5::PredType::NATIVE_DOUBLE);
//    mtype1.insertMember( "positionY", HOFFSET(s1_t, positionY), H5::PredType::NATIVE_DOUBLE);
//    mtype1.insertMember( "positionZ", HOFFSET(s1_t, positionZ), H5::PredType::NATIVE_DOUBLE);
//    mtype1.insertMember( "velocityX", HOFFSET(s1_t, velocityX), H5::PredType::NATIVE_DOUBLE);
//    mtype1.insertMember( "velocityY", HOFFSET(s1_t, velocityY), H5::PredType::NATIVE_DOUBLE);
//    mtype1.insertMember( "velocityZ", HOFFSET(s1_t, velocityZ), H5::PredType::NATIVE_DOUBLE);
//    mtype1.insertMember( "forceX", HOFFSET(s1_t, forceX), H5::PredType::NATIVE_DOUBLE);
//    mtype1.insertMember( "forceY", HOFFSET(s1_t, forceY), H5::PredType::NATIVE_DOUBLE);
//    mtype1.insertMember( "forceZ", HOFFSET(s1_t, forceZ), H5::PredType::NATIVE_DOUBLE);
//    mtype1.insertMember( "cellID", HOFFSET(s1_t, cellID), H5::PredType::NATIVE_INT);

//    /*
//      * Create the dataset.
//      */
//    H5::DataSet* dataset;
//    dataset = new H5::DataSet(file->createDataSet("ArrayOfStructures", mtype1, space));

//    /*
//      * Write data to the dataset;
//      */
//    dataset->write( s1, mtype1 );

//    /*
//      * Release resources
//      */
//    delete dataset;
//    delete file;
//    return true;
//}

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
