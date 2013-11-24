#ifndef VELOCITYVERLETINTEGRATOR_H
#define VELOCITYVERLETINTEGRATOR_H

#include <integrator/integrator.h>

class VelocityVerletIntegrator : public Integrator
{
public:
    VelocityVerletIntegrator(MoleculeSystem *moleculeSystem);
    void initialize();
    void stepForward();
};

#endif // VELOCITYVERLETINTEGRATOR_H
