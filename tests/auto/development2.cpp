#include <unittest++/UnitTest++.h>

#include "moleculesystem.h"
#include "moleculesystemcell.h"
#include "generator.h"
#include "atomtype.h"
#include "force/lennardjonesforce.h"
#include "force/fannthreeparticleforce.h"
#include "force/fanntwoparticleforce.h"
#include "modifier/andersenthermostat.h"
#include "modifier/berendsenthermostat.h"
#include "modifier/friction.h"
#include "atom.h"
#include "processor.h"
#include "utils/logging.h"
#include "force/threeparticleforce.h"
#include "integrator/velocityverletintegrator.h"
#include "filemanager.h"
#include <boost/mpi.hpp>

SUITE(Development2) {
    TEST(Dummy) {

    }
    TEST(FannForceHydrogen)
    {
        LOG(INFO) << "FannForceHydrogen test started";


        double electronMassPerAtomicMass = 1822.8896;
        // Particle types
        AtomType hydrogenType(0);
        hydrogenType.setMass(1.008*electronMassPerAtomicMass);
        hydrogenType.setNumber(1);
        hydrogenType.setAbbreviation("H");

        vector<AtomType> particleTypes;
        particleTypes.push_back(hydrogenType);

        Atom *hydrogenAtom1 = new Atom(hydrogenType);
        hydrogenAtom1->setPosition(Vector3(0.0, 0.0, 0.0));
        hydrogenAtom1->setID(1);
        Atom *hydrogenAtom2 = new Atom(hydrogenType);
        hydrogenAtom2->setPosition(Vector3(0.0, 2.0, 0.0));
        hydrogenAtom2->setID(2);

        double cutoffRadius = 10.0;

        FannTwoParticleForce testForce2;
        testForce2.setNewtonsThirdLawEnabled(true);
        testForce2.setCutoffRadius(cutoffRadius);
//        testForce2.addNetwork(hydrogenType, hydrogenType,
//                              "/home/svenni/Dropbox/studies/master/results/fann_train/20140427-141531/fann_network.net",
//                              "/home/svenni/Dropbox/studies/master/results/fann_train/20140427-141531/bounds.fann");
        testForce2.addNetwork(hydrogenType, hydrogenType,
                              "/home/svenni/Dropbox/studies/master/results/fann_train/20140507-204601/fann_network.net",
                              "/home/svenni/Dropbox/studies/master/results/fann_train/20140507-204601/bounds.fann");
        testForce2.addNetwork(hydrogenType, hydrogenType,
                              "/home/svenni/Dropbox/studies/master/results/fann_train/20140507-200839/fann_network.net",
                              "/home/svenni/Dropbox/studies/master/results/fann_train/20140507-200839/bounds.fann");

        ofstream outFile("/tmp/forcedata.out");

        rowvec xValues = linspace(0.5, 20, 1000);
        for(double x : xValues) {
            hydrogenAtom1->clearForcePotentialPressure();
            hydrogenAtom2->clearForcePotentialPressure();
            hydrogenAtom2->setPosition(Vector3(x, 0, 0));
            testForce2.calculateAndApplyForce(hydrogenAtom1, hydrogenAtom2);

            outFile << x << " " << hydrogenAtom1->potential() * 2 << " " << hydrogenAtom1->force().x() << endl;
        }
    }
    TEST(FannForceHydrogen3)
    {
        LOG(INFO) << "FannForceHydrogen test started";


        double electronMassPerAtomicMass = 1822.8896;
        // Particle types
        AtomType hydrogenType(0);
        hydrogenType.setMass(1.008*electronMassPerAtomicMass);
        hydrogenType.setNumber(1);
        hydrogenType.setAbbreviation("H");

        vector<AtomType> particleTypes;
        particleTypes.push_back(hydrogenType);

        Atom *hydrogenAtom1 = new Atom(hydrogenType);
        hydrogenAtom1->setPosition(Vector3(0.0, 0.0, 0.0));
        hydrogenAtom1->setID(1);
        Atom *hydrogenAtom2 = new Atom(hydrogenType);
        hydrogenAtom2->setPosition(Vector3(0.0, 1.4, 0.0));
        hydrogenAtom2->setID(2);
        Atom *hydrogenAtom3 = new Atom(hydrogenType);
        hydrogenAtom3->setPosition(Vector3(1.0, 0.0, 0.0));
        hydrogenAtom3->setID(3);

        double cutoffRadius = 10.0;

        FannTwoParticleForce testForce2;
        testForce2.setNewtonsThirdLawEnabled(true);
        testForce2.setCutoffRadius(cutoffRadius);
//        testForce2.addNetwork(hydrogenType, hydrogenType,
//                              "/home/svenni/Dropbox/studies/master/results/fann_train/20140427-141531/fann_network.net",
//                              "/home/svenni/Dropbox/studies/master/results/fann_train/20140427-141531/bounds.fann");
        testForce2.addNetwork(hydrogenType, hydrogenType,
                              "/home/svenni/Dropbox/studies/master/results/fann_train/20140507-204601/fann_network.net",
                              "/home/svenni/Dropbox/studies/master/results/fann_train/20140507-204601/bounds.fann");
        testForce2.addNetwork(hydrogenType, hydrogenType,
                              "/home/svenni/Dropbox/studies/master/results/fann_train/20140507-200839/fann_network.net",
                              "/home/svenni/Dropbox/studies/master/results/fann_train/20140507-200839/bounds.fann");

        FannThreeParticleForce testForce3;
        testForce3.setNewtonsThirdLawEnabled(true);
        testForce3.setCutoffRadius(cutoffRadius);
        testForce3.loadNetwork("/home/svenni/Dropbox/studies/master/results/fann_train/20140512-200948/fann_network.net",
                               "/home/svenni/Dropbox/studies/master/results/fann_train/20140512-200948/bounds.fann",
                               5.0);
        testForce3.loadNetwork("/home/svenni/Dropbox/studies/master/results/fann_train/20140512-182447/fann_network.net",
                               "/home/svenni/Dropbox/studies/master/results/fann_train/20140512-182447/bounds.fann");

        ofstream outFile("/tmp/forcedata3.out");

        rowvec xValues = linspace(0.5, 20, 1000);
        for(double x : xValues) {
            hydrogenAtom1->clearForcePotentialPressure();
            hydrogenAtom2->clearForcePotentialPressure();
            hydrogenAtom3->clearForcePotentialPressure();
            hydrogenAtom3->setPosition(Vector3(1.4 + x, 0, 0));
            testForce3.calculateAndApplyForce(hydrogenAtom1, hydrogenAtom2, hydrogenAtom3);
            testForce3.calculateAndApplyForce(hydrogenAtom2, hydrogenAtom1, hydrogenAtom3);
            testForce3.calculateAndApplyForce(hydrogenAtom3, hydrogenAtom1, hydrogenAtom2);

            outFile << x << " " << hydrogenAtom1->potential() * 2 << " " << hydrogenAtom3->force().x() << endl;
        }
    }
}
