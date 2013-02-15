#include "filemanager.h"

#include <src/moleculesystem.h>
#include <src/moleculesystemcell.h>
#include <src/atom.h>
#include <src/molecule.h>
#include <src/atomtype.h>
#include <src/integrator/integrator.h>

// System includes
#include <iomanip>
#include <armadillo>
#include <H5Cpp.h>
#include <H5File.h>
//#include <hdf5.h>

using namespace std;
using namespace arma;
using namespace H5;

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
    m_unitMass(1)
{
}

void FileManager::load(string fileName) {
    cerr << "Cannot load " << fileName << " because loading is not yet implemented!" << endl;
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
    double time = step * m_moleculeSystem->integrator()->timeStep();
    double timeStep = m_moleculeSystem->integrator()->timeStep();
    int nAtoms = m_moleculeSystem->atoms().size();
    double temperature = m_moleculeSystem->temperature();
    rowvec averageDisplacement = m_moleculeSystem->averageDisplacement();
    double averageSquareDisplacement = m_moleculeSystem->averageSquareDisplacement();
    outFile.write((char*)&step, sizeof(int));
    outFile.write((char*)&time, sizeof(double));
    outFile.write((char*)&timeStep, sizeof(double));
    outFile.write((char*)&nAtoms, sizeof(int));
    outFile.write((char*)&temperature, sizeof(double));
    outFile.write((char*)&averageDisplacement(0), sizeof(double));
    outFile.write((char*)&averageDisplacement(1), sizeof(double));
    outFile.write((char*)&averageDisplacement(2), sizeof(double));
    outFile.write((char*)&averageSquareDisplacement, sizeof(double));

    // Write atom data
    for(MoleculeSystemCell* cell : m_moleculeSystem->cells()) {
        for(Molecule* molecule : cell->molecules()) {
            for(Atom* atom : molecule->atoms()) {
                rowvec position = atom->position() * m_unitLength;
                rowvec velocity = atom->velocity() * (m_unitLength / m_unitTime);
                rowvec force = atom->force() * (m_unitMass * m_unitLength / (m_unitTime * m_unitTime));
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
    DataSpace space( 1, dim );

    /*
      * Create the file.
      */
    H5File* file = new H5File( outFileNameLocal, H5F_ACC_TRUNC );

    /*
      * Create the memory datatype.
      */
    H5::StrType string_type(H5::PredType::C_S1, 3);
    CompType mtype1( sizeof(s1_t) );
    mtype1.insertMember( "atomType", HOFFSET(s1_t, atomType), string_type);
    mtype1.insertMember( "positionX", HOFFSET(s1_t, positionX), PredType::NATIVE_DOUBLE);
    mtype1.insertMember( "positionY", HOFFSET(s1_t, positionY), PredType::NATIVE_DOUBLE);
    mtype1.insertMember( "positionZ", HOFFSET(s1_t, positionZ), PredType::NATIVE_DOUBLE);
    mtype1.insertMember( "velocityX", HOFFSET(s1_t, velocityX), PredType::NATIVE_DOUBLE);
    mtype1.insertMember( "velocityY", HOFFSET(s1_t, velocityY), PredType::NATIVE_DOUBLE);
    mtype1.insertMember( "velocityZ", HOFFSET(s1_t, velocityZ), PredType::NATIVE_DOUBLE);
    mtype1.insertMember( "forceX", HOFFSET(s1_t, forceX), PredType::NATIVE_DOUBLE);
    mtype1.insertMember( "forceY", HOFFSET(s1_t, forceY), PredType::NATIVE_DOUBLE);
    mtype1.insertMember( "forceZ", HOFFSET(s1_t, forceZ), PredType::NATIVE_DOUBLE);
    mtype1.insertMember( "cellID", HOFFSET(s1_t, cellID), PredType::NATIVE_INT);

    /*
      * Create the dataset.
      */
    DataSet* dataset;
    dataset = new DataSet(file->createDataSet("ArrayOfStructures", mtype1, space));

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

void FileManager::setOutFileName(string fileName)
{
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
