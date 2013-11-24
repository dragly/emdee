#include <generator.h>
#include <atom.h>

#include <unittest++/UnitTest++.h>

/*!
 * \brief BoltzmannDistributeVelocitiesTest ensures that the velocity distribution has a total velocity drift of zero.
 */
//TEST(BoltzmannDistributeVelocitiesTest)
//{
//    Generator generator;

//    vector<Atom*> atoms = generator.generateFcc(2,8,AtomType::argon());
//    generator.boltzmannDistributeVelocities(100, atoms);

//    Vector3 totalVelocity; // = zeros<Vector3>(3);
//    for(Atom* atom : atoms) {
//        totalVelocity += atom->velocity();
//    }
//    double maxVelocity = 0;
//    CHECK(fabs(maxVelocity) < 1e-5);
//}
