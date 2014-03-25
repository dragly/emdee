#include <iostream>
#include <sstream>
#include <fstream>
#include <fann.h>
#include <mpi/mpi.h>

using namespace std;

void train(int rank) {
    const unsigned int num_input = 3;
    const unsigned int num_output = 1;
    const unsigned int num_layers = 10;
    const unsigned int num_neurons_hidden = 12;
    const double desired_error = (const double) 0.0000000005;
    const unsigned int max_epochs = 30000;
    const unsigned int max_neurons = 1;
    const unsigned int neurons_between_reports = 1;
    const unsigned int epochs_between_reports = 1000;
    double bestTestResult = 999999999;

    bool cascade = true;
    struct fann *ann;
    if(cascade) {
        ann = fann_create_shortcut(num_layers,
                                   num_input,
                                   num_neurons_hidden, num_neurons_hidden, num_neurons_hidden,num_neurons_hidden,num_neurons_hidden,num_neurons_hidden,num_neurons_hidden,num_neurons_hidden,
                                   num_output);
    } else {
        ann = fann_create_standard(num_layers,
                                   num_input,
                                   num_neurons_hidden, num_neurons_hidden, num_neurons_hidden,num_neurons_hidden,num_neurons_hidden,num_neurons_hidden,num_neurons_hidden,num_neurons_hidden,
                                   num_output);
    }
    //    fann_set_activation_function_hidden(ann, FANN_SIGMOID_SYMMETRIC);
    //    fann_set_activation_function_output(ann, FANN_LINEAR);
    fann_set_training_algorithm(ann, FANN_TRAIN_RPROP);
    stringstream inFileName;
    inFileName << "/home/svenni/Dropbox/studies/master/code/hartree-fock/hartree-fock/tools/train.fann";
    struct fann_train_data *trainData = fann_read_train_from_file(inFileName.str().c_str());
    stringstream inFileTestName;
    inFileTestName << "/home/svenni/Dropbox/studies/master/code/hartree-fock/hartree-fock/tools/test.fann";
    struct fann_train_data *testData = fann_read_train_from_file(inFileTestName.str().c_str());
//    enum fann_activationfunc_enum activation = FANN_SIGMOID_SYMMETRIC;
//    fann_set_cascade_activation_functions(ann, &activation, 1);
    //        fann_set_cascade_num_candidate_groups(ann, 32);
    //        fann_set_cascade_candidate_limit(ann, 1000);
    //        fann_set_cascade_output_change_fraction(ann, 0.2);
//    fann_set_cascade_weight_multiplier(ann, 100.0);

    for(int i = 1; i < 20; i++) {
        cout << "---------- RANK " << rank << " ITERATION " << i << "-----------" << endl;
        //            fann_set_cascade_candidate_limit(ann, i*1000);
        fann_cascadetrain_on_data(ann, trainData, 1, neurons_between_reports, desired_error);
//        cout << "Reset and test" << endl;
        fann_reset_MSE(ann);
        fann_test_data(ann, testData);
        double testResult = fann_get_MSE(ann);
        cout << "Results from testing data: " << testResult << endl;
        if(testResult / bestTestResult > 1.001) {
            cout << "No longer converging. Stopping to avoid overfitting." << endl;
            break;
        } else {
            bestTestResult = testResult;
            stringstream outFileName;
            outFileName << "fann_network" << rank << ".net";
            fann_save(ann, outFileName.str().c_str());
            stringstream outFileName2;
            outFileName2 << "testresult_" << rank << ".out";
            ofstream testResultFile(outFileName2.str());
            testResultFile << bestTestResult;
            testResultFile.close();
        }
    }
    //    }

    fann_destroy_train(trainData);
    fann_destroy_train(testData);
    fann_destroy(ann);
}

int main(int argc, char* argv[])
{
    int rank;
    int nProcessors;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nProcessors);
    int padder = 1000;
    for(int i = 0; i < 10; i++) {
        cout << "--------------- RANK " << rank << " SET " << i << " -----------------" << endl;
        int id = rank * padder + i;
        train(id);
    }
    MPI_Finalize();
    return 0;
}

