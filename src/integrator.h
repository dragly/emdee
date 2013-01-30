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
private:
    MoleculeSystem *m_moleculeSystem;

    double dt;
};

#endif // INTEGRATOR_H
