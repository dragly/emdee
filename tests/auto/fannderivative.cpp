#include <unittest++/UnitTest++.h>

#include "moleculesystem.h"
#include "moleculesystemcell.h"
#include "generator.h"
#include "atomtype.h"
#include "force/lennardjonesforce.h"
#include "force/fannthreeparticleforce.h"
#include "force/fanntwoparticleforce.h"
#include "modifier/andersenthermostat.h"
#include "modifier/berendsenthermostat.h"
#include "modifier/friction.h"
#include "atom.h"
#include "processor.h"
#include "utils/logging.h"
#include "force/threeparticleforce.h"
#include "integrator/velocityverletintegrator.h"
#include "filemanager.h"
#include "utils/fannderivative.h"

#include <iomanip>
#include <doublefann.h>

using namespace std;

SUITE(FannDerivative) {
    TEST(Dummy) {

    }
    TEST(FannDerivative)
    {
        cout << setprecision(20);
//        string inFileName = "/home/svenni/Dropbox/studies/master/results/fann_train/20140418-165800/fann_network.net";
//        struct fann* ann = fann_create_from_file(inFileName.c_str());

        uint layerCount = 4;
        uint inputCount = 3;
        uint outputCount = 3;

        struct fann* ann = fann_create_standard(layerCount, inputCount, 5, 5, outputCount);
        int layer_number = 1;
        for(fann_layer* layer_it = ann->first_layer + 1; layer_it < ann->last_layer; layer_it++) {
            int neuron_number = 0;
            for(fann_neuron* neuron_it = layer_it->first_neuron; neuron_it < layer_it->last_neuron; neuron_it++) {
                fann_activationfunc_enum activation_function = fann_get_activation_function(ann, layer_number, neuron_number);
                switch (activation_function)
                {
                case FANN_SIGMOID_STEPWISE:
                    fann_set_activation_function(ann, FANN_SIGMOID, layer_number, neuron_number);
                    break;
                case FANN_SIGMOID_SYMMETRIC_STEPWISE:
                    fann_set_activation_function(ann, FANN_SIGMOID_SYMMETRIC, layer_number, neuron_number);
                    break;
                case FANN_THRESHOLD:
                    fann_error(NULL, FANN_E_CANT_TRAIN_ACTIVATION);
                    break;
                case FANN_THRESHOLD_SYMMETRIC:
                    fann_error(NULL, FANN_E_CANT_TRAIN_ACTIVATION);
                    break;
                default:
                    break;
                }
                neuron_number++;
            }
            layer_number++;
        }

        fann_type input[inputCount];
        for(int i = 0; i < int(inputCount); i++) {
            input[i] = 0.0;
        }

        int testCount = 100;
        for(int i = 0; i < testCount; i++) {
            fann_type value = i * 0.8 / testCount + 0.1;
            input[0] = value;
            // Numerical derivative
            fann_type h = 1e-5;
            input[0] = value + h;
            fann_type outputPlus = fann_run(ann, input)[0];
            input[0] = value - h;
            fann_type outputMinus = fann_run(ann, input)[0];
            fann_type numericDerivative = (outputPlus - outputMinus) / (2*h);
            input[0] = value;

            fann_run(ann, input);
            FannDerivative::backpropagateDerivative(ann, 0);
            CHECK_CLOSE(numericDerivative, ann->train_errors[0], 1e-4);
        }


        fann_destroy(ann);
    }
}
