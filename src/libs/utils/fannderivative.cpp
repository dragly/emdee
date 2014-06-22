#include "fannderivative.h"
#include <iostream>
#include <string.h>
#include <doublefann.h>

using namespace std;

/*!
 * \class FannDerivative
 * \brief The FannDerivative class is a helper class that is used to calculate the
 * gradient of the FANN network output with respect to its input.
 */


FannDerivative::FannDerivative()
{
}

double FannDerivative::activationDerived(unsigned int activation_function,
                                         fann_type steepness, fann_type value, fann_type sum)
{
    switch (activation_function)
    {
    case FANN_LINEAR:
    case FANN_LINEAR_PIECE:
    case FANN_LINEAR_PIECE_SYMMETRIC:
        return (fann_type) fann_linear_derive(steepness, value);
    case FANN_SIGMOID:
        return (fann_type) fann_sigmoid_derive(steepness, value);
    case FANN_SIGMOID_STEPWISE:
        fann_error(NULL, FANN_E_CANT_TRAIN_ACTIVATION);
    case FANN_SIGMOID_SYMMETRIC:
        return (fann_type) fann_sigmoid_symmetric_derive(steepness, value);
    case FANN_SIGMOID_SYMMETRIC_STEPWISE:
        fann_error(NULL, FANN_E_CANT_TRAIN_ACTIVATION);
    case FANN_GAUSSIAN:
        return (fann_type) fann_gaussian_derive(steepness, value, sum/steepness);
    case FANN_GAUSSIAN_SYMMETRIC:
        return (fann_type) fann_gaussian_symmetric_derive(steepness, value, sum/steepness);
    case FANN_ELLIOT:
        return (fann_type) fann_elliot_derive(steepness, value, sum);
    case FANN_ELLIOT_SYMMETRIC:
        return (fann_type) fann_elliot_symmetric_derive(steepness, value, sum);
    case FANN_SIN_SYMMETRIC:
        return (fann_type) fann_sin_symmetric_derive(steepness, sum/steepness);
    case FANN_COS_SYMMETRIC:
        return (fann_type) fann_cos_symmetric_derive(steepness, sum/steepness);
    case FANN_SIN:
        return (fann_type) fann_sin_derive(steepness, sum/steepness);
    case FANN_COS:
        return (fann_type) fann_cos_derive(steepness, sum/steepness);
    case FANN_THRESHOLD:
        fann_error(NULL, FANN_E_CANT_TRAIN_ACTIVATION);
    }
    return 0;
}

void FannDerivative::backpropagateDerivative(struct fann *ann, uint outputIndex)
{
    fann_type tmp_error;
    unsigned int i;
    struct fann_layer *layer_it;
    struct fann_neuron *neuron_it, *last_neuron;
    struct fann_neuron **connections;

    fann_type *error_prev_layer;
    fann_type *weights;
    const struct fann_neuron *first_neuron = ann->first_layer->first_neuron;
    struct fann_layer *last_layer = ann->last_layer;

    /* if no room allocated for the error variabels, allocate it now */
    if(ann->train_errors == NULL)
    {
        ann->train_errors = (fann_type *) calloc(ann->total_neurons, sizeof(fann_type));
        if(ann->train_errors == NULL)
        {
            fann_error((struct fann_error *) ann, FANN_E_CANT_ALLOCATE_MEM);
            return;
        }
    }
    /* clear the error variabels */
    memset(ann->train_errors, 0, (ann->total_neurons) * sizeof(fann_type));
    fann_type *error_begin = ann->train_errors;

    /* Set the error variable of the output neuron */
    fann_neuron *output_neuron = (ann->last_layer-1)->first_neuron + outputIndex;
    ann->train_errors[output_neuron - first_neuron] = activationDerived(output_neuron->activation_function,
                                                                        output_neuron->activation_steepness,
                                                                        output_neuron->value,
                                                                        output_neuron->sum);

    /* go through all the layers, from last to first.
     * And propagate the error backwards */
    for(layer_it = last_layer - 1; layer_it > ann->first_layer; --layer_it)
    {
        last_neuron = layer_it->last_neuron;

        /* for each connection in this layer, propagate the error backwards */
        if(ann->connection_rate >= 1)
        {
            if(ann->network_type == FANN_NETTYPE_LAYER)
            {
                error_prev_layer = error_begin + ((layer_it - 1)->first_neuron - first_neuron);
            }
            else
            {
                error_prev_layer = error_begin;
            }

            for(neuron_it = layer_it->first_neuron; neuron_it != last_neuron; neuron_it++)
            {
                tmp_error = error_begin[neuron_it - first_neuron];
                weights = ann->weights + neuron_it->first_con;
                for(i = neuron_it->last_con - neuron_it->first_con; i--;)
                {
                    /*printf("i = %d\n", i);
                     * printf("error_prev_layer[%d] = %f\n", i, error_prev_layer[i]);
                     * printf("weights[%d] = %f\n", i, weights[i]); */
                    error_prev_layer[i] += tmp_error * weights[i];
                }
            }
        }
        else
        {
            for(neuron_it = layer_it->first_neuron; neuron_it != last_neuron; neuron_it++)
            {

                tmp_error = error_begin[neuron_it - first_neuron];
                weights = ann->weights + neuron_it->first_con;
                connections = ann->connections + neuron_it->first_con;
                for(i = neuron_it->last_con - neuron_it->first_con; i--;)
                {
                    error_begin[connections[i] - first_neuron] += tmp_error * weights[i];
                }
            }
        }

        /* then calculate the actual errors in the previous layer */
        error_prev_layer = error_begin + ((layer_it - 1)->first_neuron - first_neuron);
        last_neuron = (layer_it - 1)->last_neuron;

        if(layer_it - 1 > ann->first_layer) {
            for(neuron_it = (layer_it - 1)->first_neuron; neuron_it != last_neuron; neuron_it++)
            {
                *error_prev_layer *= activationDerived(neuron_it->activation_function,
                                                       neuron_it->activation_steepness,
                                                       neuron_it->value,
                                                       neuron_it->sum);
                error_prev_layer++;
            }
        }
    }
}
