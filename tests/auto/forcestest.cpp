#include <unittest++/UnitTest++.h>
#include <src/atom.h>
#include <src/math/vector3.h>
#include <src/force/lennardjonesforce.h>
#include <src/force/vashishtatwoparticleforce.h>
#include <src/force/vashishtathreeparticleforce.h>

TEST(ArgonForceTest) {
    Atom* atom1 = new Atom(AtomType::argon());
    Atom* atom2 = new Atom(AtomType::argon());
    Vector3 position1(-1,0,0);
    Vector3 position2(1,0,0);
    atom1->setPosition(position1);
    atom2->setPosition(position2);
    LennardJonesForce force;
    force.setPotentialConstant(3);
    force.setNewtonsThirdLawEnabled(true);
    force.calculateAndApplyForce(atom1, atom2);
    CHECK_EQUAL(atom1->force(), -atom2->force());
    Vector3 newForce(-2977.22, 0, 0);
    CHECK_ARRAY_CLOSE(atom1->force(), newForce, 3, 0.005);
}

TEST(VashishtaForceTest) {
    AtomType silicon;
    silicon.setMass(28.0855);
    silicon.setId(14);
    silicon.setEffectiveCharge(1.6);
    silicon.setElectronicPolarizability(0.0);
    AtomType oxygen;
    oxygen.setMass(15.9994);
    oxygen.setId(8);
    oxygen.setEffectiveCharge(-0.8);
    oxygen.setElectronicPolarizability(2.40);


    Atom atom1; // = new Atom(silicon);
    Atom atom2; // = new Atom(oxygen);
    uint nSteps = 10000;
    vec positionX = linspace(100,1.15,nSteps);
    vec forceX = zeros(nSteps);
    vec potential = zeros(nSteps);
    Vector3 position1(0,0,0);
    VashishtaTwoParticleForce force;
    force.setNewtonsThirdLawEnabled(true);
    mat toPrint = zeros(nSteps, 7);
    toPrint.col(0) = positionX;
    for(int i = 0; i < 3; i++) {
        switch(i) {
        case 0:
            atom1 = Atom(silicon);
            atom2 = Atom(oxygen);
            break;
        case 1:
            atom1 = Atom(silicon);
            atom2 = Atom(silicon);
            break;
        case 2:
            atom1 = Atom(oxygen);
            atom2 = Atom(oxygen);
            break;
        default:
            throw exception();
            break;
        }
        for(uint i = 0; i < positionX.size(); i++) {
            atom1.clearForcePotentialPressure();
            atom2.clearForcePotentialPressure();
            Vector3 position2(positionX[i],0,0);
            atom1.setPosition(position1);
            atom2.setPosition(position2);
        //    force.setPotentialConstant(3);
            try {
                force.calculateAndApplyForce(&atom1, &atom2);
            } catch(exception) {

            }

            forceX[i] = atom2.force()[0];
            potential[i] = atom2.potential() + atom1.potential();
//            CHECK_EQUAL(atom2.force(), -atom1.force());
//            Vector3 newForce(-2977.22, 0, 0);
        }
        toPrint.col(i + 1) = forceX;
        toPrint.col(i + 4) = potential;
    }
    toPrint.save("output.dat", raw_ascii);
//    CHECK_ARRAY_CLOSE(atom1.force(), newForce, 3, 0.005);
}

TEST(VashishtaThreeForceTest) {
    AtomType silicon;
    silicon.setMass(28.0855);
    silicon.setId(14);
    silicon.setEffectiveCharge(1.6);
    silicon.setElectronicPolarizability(0.0);
    AtomType oxygen;
    oxygen.setMass(15.9994);
    oxygen.setId(8);
    oxygen.setEffectiveCharge(-0.8);
    oxygen.setElectronicPolarizability(2.40);


    Atom atom1; // = new Atom(silicon);
    Atom atom2; // = new Atom(oxygen);
    Atom atom3; // = new Atom(oxygen);
    uint nSteps = 10000;
    vec positionX = linspace(100,1.15,nSteps);
    vec forceX = zeros(nSteps);
    vec potential = zeros(nSteps);
    Vector3 position1(0,0,0);
    VashishtaThreeParticleForce force;
    force.setNewtonsThirdLawEnabled(true);

    atom1.setPosition(Vector3(0,0,0));
    atom2.setPosition(Vector3(2,0,0));
    atom3.setPosition(Vector3(0,3,0));

    force.calculateAndApplyForce(&atom1, &atom2, &atom3);

//    CHECK_ARRAY_CLOSE(atom1.force(), newForce, 3, 0.005);
}
