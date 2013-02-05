#include "generator.h"

#include "atom.h"
#include "molecule.h"
/*!
 * \brief Generator::Generator has multiple functions to generate different setups of atoms.
 * Something is longer than shorter is simple. This is bla bla bla boom.
 */

Generator::Generator() :
    m_unitLength(1)
{
}

void Generator::loadConfiguration(Config *config)
{
    m_unitLength = config->lookup("units.length");
}

vector<Molecule*> Generator::generateFcc(double b, int nCells, AtomType atomType) {
    vector<Molecule*> moleculeList;
    rowvec offset = zeros<rowvec>(3);
    double bUnit = b; // / m_unitLength;
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
                        face << bUnit / 2 << bUnit / 2 <<  0;
                        break;
                    case 2:
                        face << bUnit / 2 << 0 << bUnit / 2;
                        break;
                    case 3:
                        face << 0 << bUnit / 2 << bUnit / 2;
                        break;
                    }
                    rowvec position = zeros<rowvec>(3);
                    position = offset + face;
                    molecule->setPosition(position);
                    moleculeList.push_back(molecule);
                }
                offset(2) += bUnit;
            }
            offset(1) += bUnit;
            offset(2) = 0;
        }
        offset(0) += bUnit;
        offset(1) = 0;
    }
    m_lastBoundaries << 0 << 0 << 0 << endr
                        << nCells * bUnit << nCells * bUnit << nCells * bUnit;

    cout << "Generated " << moleculeList.size() << " molecules in FCC structure!" << endl;
    return moleculeList;
}

void Generator::boltzmannDistributeVelocities(vector<Molecule*> molecules) {
    for(Molecule* molecule : molecules) {
        rowvec velocity = 0.1 * randn<rowvec>(3); // / m_unitLength;
        molecule->setVelocity(velocity);
    }
    cout << "Boltzmann distributed velocities for " << molecules.size() << " molecules!" << endl;
}

void Generator::setUnitLength(double unitLength)
{
    m_unitLength = unitLength;
}

const mat &Generator::lastBoundaries() const
{
    return m_lastBoundaries;
}
