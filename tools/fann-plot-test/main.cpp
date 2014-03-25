#include <iostream>
#include <iomanip>
#include <fstream>
#include <armadillo>
#include <fann.h>

using namespace std;
using arma::vec;
using arma::linspace;


double rescale(double value, double valueMin, double valueMax) {
    return (value - valueMin) / (valueMax - valueMin) * 0.8 + 0.1;
}

int main()
{
    string inFileName = "/home/svenni/Dropbox/projects/programming/build-emdee-Desktop_Qt_5_2_0_with_GDB-Release/tools/fann-trainer/xor_float_3.net";
    vec angles = linspace(M_PI/3, M_PI, 10);
    vec r12s = linspace(1.0, 6.0, 10);
    vec r13s = linspace(1.0, 6.0, 10);
    struct fann *ann = fann_create_from_file(inFileName.c_str());

    for(int i = 0; i < int(r12s.n_elem); i++) {
        double r12 = r12s[i];
        stringstream outFileName;
        outFileName << "testdata" << setprecision(4) << r12;
        ofstream outFile(outFileName.str());
        for (int j = 0; j < int(r13s.n_elem); ++j) {
            double r13 = r13s[j];
            for (int k = 0; k < int(angles.n_elem); ++k) {
                double angle = angles[k];
                fann_type input[3];
                input[0] = rescale(r12, r12s.min(), r12s.max());
                input[1] = rescale(r13, r13s.min(), r13s.max());
                input[2] = rescale(angle, angles.min(), angles.max());
                fann_type *energy = fann_run(ann, input);
                outFile << energy[0] << " ";
            }
            outFile << "\n";
        }
        outFile.close();
    }
    fann_destroy(ann);
    return 0;
}

