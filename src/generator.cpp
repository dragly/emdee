#include "generator.h"

#include "atom.h"
#include "molecule.h"
/*!
 * \brief Generator::Generator has multiple functions to generate different setups of atoms.
 * Something is longer than shorter is simple. This is bla bla bla boom.
 */

Generator::Generator()
{
}

vector<Molecule*> Generator::generateFcc(int nCells, double b, AtomType atomType) {
    vector<Molecule*> moleculeList;
    rowvec offset = zeros<rowvec>(3);
    for(int i = 0; i < nCells; i++) {
        for(int j = 0; j < nCells; j++) {
            for(int k = 0; k < nCells; k++) {
                for(int atomi = 0; atomi < 4; atomi++) {
                    Molecule* molecule = new Molecule();
                    Atom* atom = new Atom(molecule, atomType);
                    molecule->addAtom(atom);
                    rowvec face = zeros<rowvec>(3);
                    switch(atomi) {
                    case 0:
                        break;
                    case 1:
                        face << b / 2 << b / 2 <<  0;
                        break;
                    case 2:
                        face << b / 2 << 0 << b / 2;
                        break;
                    case 3:
                        face << 0 << b / 2 << b / 2;
                        break;
                    }
                    rowvec position = zeros<rowvec>(3);
                    position = offset + face;
                    molecule->setPosition(position);
                    moleculeList.push_back(molecule);
                }
                offset(2) += b;
            }
            offset(1) += b;
            offset(2) = 0;
        }
        offset(0) += b;
        offset(1) = 0;
    }
    cout << "Generated " << moleculeList.size() << " molecules in FCC structure!" << endl;
    return moleculeList;
}

void Generator::boltzmannDistributeVelocities(vector<Molecule*> molecules) {
    for(Molecule* molecule : molecules) {
        rowvec velocity = randn<rowvec>(3);
        molecule->setVelocity(velocity);
    }
}
