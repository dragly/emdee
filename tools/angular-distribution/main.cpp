#include <iostream>

using namespace std;

#include <src/moleculesystem.h>
#include <src/atom.h>
#include <src/math/vector3.h>
#include <fstream>
#include <armadillo>
#include <iomanip>

using namespace arma;
using namespace std;

int main(int argc, char* argv[])
{
    if(argc < 7) {
        cout << "EROR! Too few arguments..." << endl;
        cout << "Usage: angular-distribution inputfile outputfile atom1number atom2number atom3number nBins cutoff" << endl;
        exit(0);
    }
    int atomTypeNumber1 = atoi(argv[3]);
    int atomTypeNumber2 = atoi(argv[4]);
    int atomTypeNumber3 = atoi(argv[5]);
    int nBins = atoi(argv[6]);
    double cutoffDistance = atof(argv[7]);
    cout << "Starting angular distribution calculator" << endl;
    cout << "cutoffDistance: " << cutoffDistance << endl;

//    int a1 = atomTypeNumber1;
//    int a2 = atomTypeNumber2;
//    int a3 = atomTypeNumber3;
//    imat combinations(6,3);
//    irowvec combination;
//    combination = {a1, a2, a3};
//    combinations.row(0) = combination;
//    combination = {a1, a3, a2};
//    combinations.row(1) = combination;
//    combination ={a2, a1, a3};
//    combinations.row(2) = combination;
//    combination ={a2, a3, a1};
//    combinations.row(3) = combination;
//    combination ={a3, a1, a2};
//    combinations.row(4) = combination;
//    combination ={a3, a2, a1};
//    combinations.row(5) = combination;

    int a1 = atomTypeNumber1;
    int a2 = atomTypeNumber2;
    int a3 = atomTypeNumber3;
    imat combinations(2,3);
    irowvec combination;
    combination = {a1, a2, a3};
    combinations.row(0) = combination;
    combination = {a3, a2, a1};
    combinations.row(1) = combination;

    MoleculeSystem system;
    system.load(argv[1]);
    cout << "Calculating angles" << endl;


//    rowvec sideLengths = system.boundaries().row(1) - system.boundaries().row(0);
//    Vector3 distance;
//    vec angles = zeros(int(system.atoms().size() * (system.atoms().size() - 1) * (system.atoms().size() - 2) / 2.0));
//    int counter = 0;
    vec binEdges = linspace(0, M_PI, nBins + 1);
    double binSize = binEdges(1) - binEdges(0);
    vec binContent = zeros(nBins);
    for(uint i = 0; i < system.atoms().size(); i++) {
//        cout << "i: " << i << endl;
        Atom* atom1 = system.atoms()[i];
        for(uint j = 0; j < system.atoms().size(); j++) {
            if(i == j) {
                continue;
            }
//            cout << "j: " << j << endl;
            Atom* atom2 = system.atoms()[j];
            Vector3 rij = atom1->position() - atom2->position();
            double lij = sqrt(dot(rij,rij));
            if(lij > cutoffDistance) {
                continue;
            }
            for(uint k = 0; k < system.atoms().size(); k++) {
                if(i == k || j == k) {
                    continue;
                }
                Atom* atom3 = system.atoms()[k];
                Vector3 rik = atom3->position() - atom2->position();
                double lik = sqrt(dot(rik,rik));
                if(lik > cutoffDistance) {
                    continue;
                }

                bool foundCombination = false;
                irowvec combination = {atom1->type().number(),
                                    atom2->type().number(),
                                    atom3->type().number()};
                for(uint iComb = 0; iComb < combinations.n_elem; iComb++) {
                    umat equality = (combinations.row(iComb) == combination);
                    if(equality.min() == 1) {
                        foundCombination = true;
                    }
                }
                if(!foundCombination) {
                    continue;
                }

                double cosangle = dot(rij,rik) / (lij*lik);

                double angle = acos(cosangle);

//                angles(counter);
//                counter += 1;
                int bin = angle / binSize;
                binContent(bin) += 1;

                //            distance = atom2->position() - atom1->position();

                //            for(int k = 0; k < 3; k++) {
                //                if (abs(distance(k)) < abs(distance(k) + sideLengths(k)) && abs(distance(k)) < abs(distance(k) - sideLengths(k))) {
                //                    continue;
                //                } else if (abs(distance(k) + sideLengths(k)) <  abs(distance(k) - sideLengths(k))) {
                //                    distance(k) = distance(k) + sideLengths(k);
                //                } else {
                //                    distance(k) = distance(k) - sideLengths(k);
                //                }
                //            }
            }
            //            distances(counter) = sqrt(dot(distance, distance));
            //            counter += 1;
        }
    }
    cout << "Writing " << sum(binContent) << " angles to file" << endl;
    ofstream distancesFile;
    distancesFile.open(argv[2], ios::binary);
    distancesFile.write((char*)&nBins, sizeof(int));
    for(uint i = 0; i < binEdges.n_elem; i++) {
        distancesFile.write((char*)&binEdges(i), sizeof(double));
    }
    for(uint i = 0; i < binContent.n_elem; i++) {
//        cout << binContent(i) << endl;
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

