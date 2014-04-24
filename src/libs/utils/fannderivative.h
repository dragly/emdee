#ifndef FANNDERIVATIVE_H
#define FANNDERIVATIVE_H

#include <doublefann.h>

class FannDerivative
{
public:
    FannDerivative();
    static double activationDerived(unsigned int activation_function, fann_type steepness, fann_type value, fann_type sum);
    static void backpropagateDerivative(fann *ann, uint outputIndex);
};

#endif // FANNDERIVATIVE_H
