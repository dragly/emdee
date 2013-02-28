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
    vector<Atom*> atomList;
    rowvec offset = zeros<rowvec>(3);
    for(int i = 0; i < nCells; i++) {
        for(int j = 0; j < nCells; j++) {
            for(int k = 0; k < nCells; k++) {
                for(int atomi = 0; atomi < 4; atomi++) {
                    Atom* atom = new Atom(atomType);
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
                    atomList.push_back(atom);
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

    cout << "Generated " << atomList.size() << " atoms in FCC structure with side length " << sideLength << endl;
    cout << "Boundaries are " << m_lastBoundaries << endl;
    return atomList;
}

/*!
 * \brief Generator::boltzmannDistributeVelocities
 * \param temperature unitless temperature
 * \param atoms a vector of atoms to apply the Boltzmann distribution of velocities on.
 *
 * \note The Boltzmann constant should have been baked into the unitless temperature.
 *
 */
void Generator::boltzmannDistributeVelocities(double temperature, const vector<Atom*>& atoms) {
    double averageVelocity = 0;
    rowvec totalVelocity = zeros<rowvec>(m_nDimensions);
    for(Atom* atom : atoms) {
        rowvec velocity = randn<rowvec>(m_nDimensions);
        velocity *= sqrt(temperature / atom->mass());
        totalVelocity += velocity;
        atom->setVelocity(velocity);
    }
    // Remove total linear momentum
    rowvec velocityToRemove = totalVelocity / atoms.size();
    averageVelocity = 0;
    for(Atom* atom : atoms) {
        rowvec newVelocity = atom->velocity() - velocityToRemove;
        atom->setVelocity(newVelocity);
        averageVelocity += norm(newVelocity, 2) / atoms.size();
        totalVelocity += newVelocity;
    }
    cout << "Boltzmann distributed velocities for " << atoms.size() << " atoms!" << endl;
    cout << "Average velocity is " << averageVelocity << endl;
}

void Generator::uniformDistributeVelocities(double maxVelocity, vector<Atom*> atoms) {

    rowvec totalVelocity = zeros<rowvec>(m_nDimensions);
    for(Atom* atom : atoms) {
        rowvec velocity = randu<rowvec>(m_nDimensions);
        velocity -= 0.5 * ones<rowvec>(m_nDimensions);
        velocity *= 2 * maxVelocity;
        totalVelocity += velocity;
        atom->setVelocity(velocity);
    }
    // Remove total linear momentum
    rowvec velocityToRemove = totalVelocity / atoms.size();
    totalVelocity.zeros();
    for(Atom* atom : atoms) {
        rowvec newVelocity = atom->velocity() - velocityToRemove;
        atom->setVelocity(newVelocity);
    }
    cout << "Uniformly distributed velocities for " << atoms.size() << " atoms!" << endl;
}

const mat &Generator::lastBoundaries() const
{
    return m_lastBoundaries;
}
