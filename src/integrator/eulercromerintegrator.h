#ifndef EULERCROMERINTEGRATOR_H
#define EULERCROMERINTEGRATOR_H

#include <src/integrator/integrator.h>

class EulerCromerIntegrator : public Integrator
{
public:
    EulerCromerIntegrator(MoleculeSystem *moleculeSystem);
    void initialize();
    void stepForward();
};

#endif // EULERCROMERINTEGRATOR_H
