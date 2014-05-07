#include "moleculesystem.h"
#include "atom.h"
#include "math/vector3.h"
#include <fstream>
#include <armadillo>
#include <iomanip>
#include <iostream>
#ifdef MD_USE_MPI
#include <boost/mpi.hpp>
namespace mpi = boost::mpi;
#endif

using namespace arma;
using namespace std;

int main(int argc, char* argv[])
{
#ifdef MD_USE_MPI
    mpi::environment env(argc, argv);
    mpi::communicator world;
    (void)env;
    (void)world;
#endif
    if(argc < 5) {
        cout << "EROR! Too few arguments..." << endl;
        cout << "Usage: radial-distribution inputfile outputfile atom1number atom2number" << endl;
        exit(0);
    }
    int atomTypeNumber1 = atoi(argv[3]);
    int atomTypeNumber2 = atoi(argv[4]);

    MoleculeSystem system;
    system.load(argv[1]);
    cout << "Calculating distances" << endl;

    int nBins = 1000;

    rowvec sideLengths = system.boundaries().row(1) - system.boundaries().row(0);
    Vector3 distance;
    vec distances = zeros(int(system.atoms().size() * (system.atoms().size() - 1) / 2.0));
    int counter = 0;
    for(uint i = 0; i < system.atoms().size(); i++) {
        Atom* atom1 = system.atoms()[i];
        for(uint j = i + 1; j < system.atoms().size(); j++) {
            Atom* atom2 = system.atoms()[j];

            if(!((atom1->type().number() == atomTypeNumber1 && atom2->type().number() == atomTypeNumber2) || (atom1->type().number() == atomTypeNumber2 && atom2->type().number() == atomTypeNumber1))) {
                continue;
            }

            distance = atom2->position() - atom1->position();

            for(int k = 0; k < 3; k++) {
                if (abs(distance(k)) < abs(distance(k) + sideLengths(k)) and abs(distance(k)) < abs(distance(k) - sideLengths(k))) {
                    continue;
                } else if (abs(distance(k) + sideLengths(k)) <  abs(distance(k) - sideLengths(k))) {
                    distance(k) = distance(k) + sideLengths(k);
                } else {
                    distance(k) = distance(k) - sideLengths(k);
                }
            }
            distances(counter) = sqrt(dot(distance, distance));
            counter += 1;
        }
    }
    vec binEdges = linspace(0, distances.max(), nBins + 1);
    double binSize = binEdges(1) - binEdges(0);
    vec binContent = zeros(nBins);
    cout << "Creating histogram " << endl;
    for(uint i = 0; i < distances.n_rows; i++) {
        double distance = distances(i);
        if(distance > 0) {
            int bin = distance / binSize;
            binContent(bin) += 1;
        }
    }
    cout << "Writing distances to file" << endl;
    ofstream distancesFile;
    distancesFile.open(argv[2], ios::binary);
    distancesFile.write((char*)&nBins, sizeof(int));
    for(uint i = 0; i < binEdges.n_elem; i++) {
        distancesFile.write((char*)&binEdges(i), sizeof(double));
    }
    for(uint i = 0; i < binContent.n_elem; i++) {
        distancesFile.write((char*)&binContent(i), sizeof(double));
    }
    distancesFile.close();
//    for(int iCol = 0; iCol < 3; iCol++) {
//        ofstream distancesFile;
//        stringstream distancesFileBase;
//        distancesFileBase << argv[2];
//        distancesFileBase << "." << setw(4) << setfill('0') << iCol;
//        distancesFile.open(distancesFileBase.str(), ios::binary);
//        for(uint i = 0; i < distances.n_rows; i++) {
//    //        cout << distances(i) << endl;
//            distancesFile.write((char*)&distances(i,iCol), sizeof(double));
//        }
//        distancesFile.close();
//    }
    return 0;
}

