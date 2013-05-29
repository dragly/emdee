#include <unittest++/UnitTest++.h>
#include <src/atom.h>
#include <src/math/vector3.h>
#include <src/force/lennardjonesforce.h>
#include <src/force/vashishtatwoparticleforce.h>
#include <src/force/vashishtathreeparticleforce.h>
#include <src/moleculesystem.h>
#include <src/filemanager.h>
#include <src/integrator/integrator.h>

void setupSixAtoms(const vector<Atom*> atoms) {
    Vector3 offset(-12, -12, -12);
    atoms[0]->setPosition(Vector3(10,12,12) + offset);
    atoms[1]->setPosition(Vector3(17.5,12.5,12.5) + offset);

    atoms[2]->setPosition(Vector3(11.5,13.6,11.5) + offset);
    atoms[3]->setPosition(Vector3(12.5,10.3,12.2) + offset);
    atoms[4]->setPosition(Vector3(15.5,13.0,12.5) + offset);
    atoms[5]->setPosition(Vector3(15,11.3,11.9) + offset);

    for(Atom* atom : atoms) {
        atom->setVelocity(Vector3(0,0,0));
        atom->clearForcePotentialPressure();
    }
}

TEST(ForceCellTest) {
    AtomType silicon(0);
    silicon.setMass(28.0855);
    silicon.setNumber(14);
    silicon.setEffectiveCharge(1.6);
    silicon.setElectronicPolarizability(0.0);
    AtomType oxygen(1);
    oxygen.setMass(15.9994);
    oxygen.setNumber(8);
    oxygen.setEffectiveCharge(-0.8);
    oxygen.setElectronicPolarizability(2.40);

    vector<Vector3> positions1;
    vector<Vector3> positions2;

    vector<Vector3> forces1;
    vector<Vector3> forces2;

    vector<Atom*> atoms;


    Atom* atom;
    for(int i = 0; i < 2; i++) {
        atom = new Atom(silicon);
        atoms.push_back(atom);
    }
    for(int i = 0; i < 4; i++) {
        atom = new Atom(oxygen);
        atoms.push_back(atom);
    }
    setupSixAtoms(atoms);

    VashishtaThreeParticleForce threeParticleForce;
    VashishtaTwoParticleForce twoParticleForce;

    MoleculeSystem system;
    system.integrator()->setTimeStep(1);
    system.addAtoms(atoms);
    system.setBoundaries(0,30);
    system.setupCells(10);
    system.setSaveEnabled(true);
    FileManager fileManager(&system);

    fileManager.setOutFileName("/tmp/data*.bin");
    system.setFileManager(&fileManager);

    system.addTwoParticleForce(&twoParticleForce);
    system.addThreeParticleForce(&threeParticleForce);

    system.setNSimulationSteps(20);
    system.simulate();
    for(Atom* atom : atoms) {
        positions1.push_back(atom->position());
        forces1.push_back(atom->force());
    }

    setupSixAtoms(atoms);
    system.setupCells(5);
    system.simulate();

    for(Atom* atom : atoms) {
        positions2.push_back(atom->position());
        forces2.push_back(atom->force());
    }

    for(uint i = 0; i < positions1.size(); i++) {
        CHECK_ARRAY_CLOSE(positions1[i], positions2[i], 3, 1e-10);
        CHECK_ARRAY_CLOSE(forces1[i], forces2[i], 3, 1e-10);
    }
}
