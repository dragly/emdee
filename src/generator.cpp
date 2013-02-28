#include <src/generator.h>

#include <src/atom.h>
#include <src/molecule.h>
/*!
 * \brief Generator::Generator has multiple functions to generate different setups of atoms.
 * Something is longer than shorter is simple. This is bla bla bla boom.
 */

Generator::Generator() :
    m_unitLength(1),
    m_nDimensions(3)
{
}

/*!
 * \brief Generator::generateFcc
 * \param b
 * \param nCells
 * \param atomType
 *
 * All molecules start out with displacement zero.
 *
 * \return
 */
vector<Atom*> Generator::generateFcc(double sideLength, int nCells, AtomType atomType) {
    vector<Atom*> moleculeList;
    rowvec offset = zeros<rowvec>(3);
    for(int i = 0; i < nCells; i++) {
        for(int j = 0; j < nCells; j++) {
            for(int k = 0; k < nCells; k++) {
                for(int atomi = 0; atomi < 4; atomi++) {
                    Atom* atom = new Atom(atomType);
//                    Atom_old* atom = new Atom_old(molecule, atomType);
//                    molecule->addAtom(atom);
                    rowvec face = zeros<rowvec>(3);
                    switch(atomi) {
                    case 0:
                        break;
                    case 1:
                        face << sideLength / 2 << sideLength / 2 <<  0;
                        break;
                    case 2:
                        face << sideLength / 2 << 0 << sideLength / 2;
                        break;
                    case 3:
                        face << 0 << sideLength / 2 << sideLength / 2;
                        break;
                    }
                    rowvec position = zeros<rowvec>(3);
                    position = offset + face;
                    atom->setPosition(position);
                    atom->clearDisplacement();
                    moleculeList.push_back(atom);
                }
                offset(2) += sideLength;
            }
            offset(1) += sideLength;
            offset(2) = 0;
        }
        offset(0) += sideLength;
        offset(1) = 0;
    }
    m_lastBoundaries << 0 << 0 << 0 << endr
                        << nCells * sideLength << nCells * sideLength << nCells * sideLength;

    cout << "Generated " << moleculeList.size() << " molecules in FCC structure with side length " << sideLength << endl;
    cout << "Boundaries are " << m_lastBoundaries << endl;
    return moleculeList;
}

/*!
 * \brief Generator::boltzmannDistributeVelocities
 * \param temperature unitless temperature
 * \param molecules a vector of molecules to apply the Boltzmann distribution of velocities on.
 *
 * \note The Boltzmann constant should have been baked into the unitless temperature.
 *
 */
void Generator::boltzmannDistributeVelocities(double temperature, const vector<Atom*>& molecules) {
    double averageVelocity = 0;
    rowvec totalVelocity = zeros<rowvec>(m_nDimensions);
    for(Atom* molecule : molecules) {
        rowvec velocity = randn<rowvec>(m_nDimensions);
        velocity *= sqrt(temperature / molecule->mass());
        totalVelocity += velocity;
        molecule->setVelocity(velocity);
    }
    // Remove total linear momentum
    rowvec velocityToRemove = totalVelocity / molecules.size();
    averageVelocity = 0;
    for(Atom* molecule : molecules) {
        rowvec newVelocity = molecule->velocity() - velocityToRemove;
        molecule->setVelocity(newVelocity);
        averageVelocity += norm(newVelocity, 2) / molecules.size();
        totalVelocity += newVelocity;
    }
    cout << "Boltzmann distributed velocities for " << molecules.size() << " molecules!" << endl;
    cout << "Average velocity is " << averageVelocity << endl;
}

void Generator::uniformDistributeVelocities(double maxVelocity, vector<Atom*> molecules) {

    rowvec totalVelocity = zeros<rowvec>(m_nDimensions);
    for(Atom* molecule : molecules) {
        rowvec velocity = randu<rowvec>(m_nDimensions);
        velocity -= 0.5 * ones<rowvec>(m_nDimensions);
        velocity *= 2 * maxVelocity;
        totalVelocity += velocity;
        molecule->setVelocity(velocity);
    }
    // Remove total linear momentum
    rowvec velocityToRemove = totalVelocity / molecules.size();
    totalVelocity.zeros();
    for(Atom* molecule : molecules) {
        rowvec newVelocity = molecule->velocity() - velocityToRemove;
        molecule->setVelocity(newVelocity);
    }
    cout << "Uniformly distributed velocities for " << molecules.size() << " molecules!" << endl;
}

const mat &Generator::lastBoundaries() const
{
    return m_lastBoundaries;
}
