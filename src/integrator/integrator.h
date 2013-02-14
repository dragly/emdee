#ifndef INTEGRATOR_H
#define INTEGRATOR_H

class MoleculeSystem;

#include <armadillo>

using namespace arma;

class Integrator
{
public:
    Integrator(MoleculeSystem *moleculeSystem);

    virtual void stepForward() = 0;
    virtual void initialize() = 0;
    void setTimeStep(double timeStep);
    double timeStep() {
        return m_timeStep;
    }

protected:
    MoleculeSystem *m_moleculeSystem;

    double m_timeStep;
};

#endif // INTEGRATOR_H
