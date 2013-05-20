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


    Atom atom1(silicon); // = new Atom(silicon);
    Atom atom2(silicon); // = new Atom(oxygen);
    Atom atom3(oxygen); // = new Atom(oxygen);
    uint nSteps = 10000;
    vec positionX = linspace(100,1.15,nSteps);
    vec forceX = zeros(nSteps);
    vec potential = zeros(nSteps);
    Vector3 position1(0,0,0);
    VashishtaThreeParticleForce force;
    force.setNewtonsThirdLawEnabled(true);

    for(int test = 0; test < 10; test++) {

        vec3 randomPosition1 = randn(3) * 0.1;
        vec3 randomPosition2 = randn(3) * 0.1;
        vec3 randomPosition3 = randn(3) * 0.1;

        atom1.setPosition(Vector3(randomPosition1));
        atom2.setPosition(Vector3(randomPosition2));
        atom3.setPosition(Vector3(randomPosition3));

        atom1.clearForcePotentialPressure();
        atom2.clearForcePotentialPressure();
        atom3.clearForcePotentialPressure();

        force.calculateAndApplyForce(&atom1, &atom2, &atom3);

//        cout << "Atom1: " << atom1.position() << "; " << atom1.force() << endl;
//        cout << "Atom2: " << atom2.position() << "; " << atom2.force() << endl;
//        cout << "Atom3: " << atom3.position() << "; " << atom3.force() << endl;

        Vector3 forceResults[3];
        forceResults[0] = atom1.force();
        forceResults[1] = atom2.force();
        forceResults[2] = atom3.force();

        for(int iAtom = 0; iAtom < 3; iAtom++) {
    //        cout << "iAtom: " << iAtom << endl;
            Atom *atom;
            switch(iAtom) {
            case 0:
                atom = &atom1;
                break;
            case 1:
                atom = &atom2;
                break;
            case 2:
                atom = &atom3;
                break;
            }

            double h = 0.00001;
            Vector3 derivativeResult;
            derivativeResult.zeros();
            Vector3 forceResult = forceResults[iAtom];
            Vector3 initialPosition = atom->position();
            for(int i = 0; i < 3; i++) {
                Vector3 minusPosition = initialPosition;
                Vector3 plusPosition = initialPosition;
                minusPosition[i] -= h;
                plusPosition[i] += h;
                atom->clearForcePotentialPressure();
                atom->setPosition(plusPosition);
                force.calculateAndApplyForce(&atom1,&atom2,&atom3);
                double plusPotential = atom->potential() * 3;
                atom->clearForcePotentialPressure();
                atom->setPosition(minusPosition);
                force.calculateAndApplyForce(&atom1,&atom2,&atom3);
                double minusPotential = atom->potential() * 3;

                double derivative = (plusPotential - minusPotential) / (2*h);
    //            cout << "derivative: " << derivative << endl;
    //            cout << "from " << "(" << plusPotential << "-" << minusPotential << ") / (2*" << h<< ")" << endl;
                derivativeResult[i] = derivative;
            }
            atom->setPosition(initialPosition);
            CHECK_ARRAY_CLOSE(derivativeResult, forceResult, 3, 0.005);
        }
    }
//    CHECK_ARRAY_CLOSE(atom1.force(), newForce, 3, 0.005);
}
