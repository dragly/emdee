#include <unittest++/UnitTest++.h>
#include <src/atom.h>
#include <src/math/vector3.h>
#include <src/force/lennardjonesforce.h>
#include <src/force/vashishtatwoparticleforce.h>

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
    Atom* atom1 = new Atom(AtomType::argon());
    Atom* atom2 = new Atom(AtomType::argon());
    Vector3 position1(-2,0,0);
    Vector3 position2(2,0,0);
    atom1->setPosition(position1);
    atom2->setPosition(position2);
    VashishtaTwoParticleForce force;
//    force.setPotentialConstant(3);
    force.setNewtonsThirdLawEnabled(true);
    force.calculateAndApplyForce(atom1, atom2);
    CHECK_EQUAL(atom1->force(), -atom2->force());
    Vector3 newForce(-2977.22, 0, 0);
//    CHECK_ARRAY_CLOSE(atom1->force(), newForce, 3, 0.005);
}
