#ifndef INTEGRATOR_H
#define INTEGRATOR_H

class MoleculeSystem;

#include <armadillo>

using namespace arma;

class Integrator
{
public:
    Integrator(MoleculeSystem *moleculeSystem);

    void stepForward();
    void initialize();
    void setTimeStep(double timeStep);
private:
    MoleculeSystem *m_moleculeSystem;

    double m_timeStep;
};

#endif // INTEGRATOR_H
