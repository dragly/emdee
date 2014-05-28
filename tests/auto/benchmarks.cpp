#include <unittest++/UnitTest++.h>

#include "moleculesystem.h"
#include "moleculesystemcell.h"
#include "generator.h"
#include "atomtype.h"
#include "force/lennardjonesforce.h"
#include "force/fannthreeparticleforce.h"
#include "force/fanntwoparticleforce.h"
#include "force/kohenthreeparticleforce.h"
#include "modifier/andersenthermostat.h"
#include "modifier/berendsenthermostat.h"
#include "modifier/friction.h"
#include "atom.h"
#include "processor.h"
#include "utils/logging.h"
#include "force/threeparticleforce.h"
#include "integrator/velocityverletintegrator.h"
#include "filemanager.h"

#ifdef MD_USE_MPI
SUITE(Benchmark) {
    TEST(Dummy) {

    }
    TEST(FannForceBenchmark) {
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
        hydrogenAtom2->setPosition(Vector3(1.4, 0.0, 0.0));
        hydrogenAtom2->setID(2);
        Atom *hydrogenAtom3 = new Atom(hydrogenType);
        hydrogenAtom3->setPosition(Vector3(0.0, 1.4, 0.0));
        hydrogenAtom3->setID(3);

        double cutoffRadius = 6.0;

        FannTwoParticleForce testForce2;
        testForce2.setNewtonsThirdLawEnabled(true);
        testForce2.setCutoffRadius(cutoffRadius);
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

        int nCalculations = 20000000;
        boost::mpi::timer timer;
        timer.restart();
        for(int i = 0; i < nCalculations; i++) {
            testForce2.calculateAndApplyForce(hydrogenAtom1, hydrogenAtom2);
        }
        cout << "Time neural network 2-body: " << timer.elapsed() << endl;

        timer.restart();
        for(int i = 0; i < nCalculations; i++) {
            testForce3.calculateAndApplyForce(hydrogenAtom1, hydrogenAtom2, hydrogenAtom3);
        }
        cout << "Time neural network 3-body: " << timer.elapsed() << endl;

        KohenThreeParticleForce kohenForce3;


        //        timer.restart();
        //        for(int i = 0; i < 1000000; i++) {
        //            testForce2.calculateAndApplyForce(hydrogenAtom1, hydrogenAtom2);
        //        }
        //        cout << hydrogenAtom1->force() << " " << hydrogenAtom2->force() << endl;
        //        cout << "Time neural network 2-body: " << timer.elapsed() << endl;


        timer.restart();
        for(int i = 0; i < nCalculations; i++) {
            kohenForce3.calculateAndApplyForce(hydrogenAtom1, hydrogenAtom2, hydrogenAtom3);
        }
        cout << "Time Kohen 3-body: " << timer.elapsed() << endl;
    }
}

#endif
