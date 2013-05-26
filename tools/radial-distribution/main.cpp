#include <iostream>

using namespace std;

#include <src/moleculesystem.h>
#include <src/atom.h>>
#include <src/math/vector3.h>
#include <fstream>
#include <armadillo>

using namespace arma;
using namespace std;

int main(int argc, char* argv[])
{
    if(argc < 3) {
        cout << "EROR! Too few arguments, needs input and output file..." << endl;
        exit(0);
    }

    MoleculeSystem system;
    system.load(argv[1]);
    cout << "Calculating distances" << endl;

    rowvec sideLengths = system.boundaries().row(1) - system.boundaries().row(0);
    Vector3 distance;
    rowvec distances = zeros<rowvec>(int(system.atoms().size() * (system.atoms().size() - 1) / 2.0));
    int counter = 0;
    for(uint i = 0; i < system.atoms().size(); i++) {
        Atom* atom1 = system.atoms()[i];
        for(uint j = i + 1; j < system.atoms().size(); j++) {
            Atom* atom2 = system.atoms()[j];
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
    cout << "Writing distances to file" << endl;
    ofstream distancesFile;
    distancesFile.open(argv[2], ios::binary);
    for(uint i = 0; i < distances.n_elem; i++) {
//        cout << distances(i) << endl;
        distancesFile.write((char*)&distances(i), sizeof(double));
    }
    distancesFile.close();
    return 0;
}

