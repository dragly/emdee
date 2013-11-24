#include <generator.h>

#include <atom.h>
/*!
 * \brief Generator::Generator has multiple functions to generate different setups of atoms.
 * Something is longer than shorter is simple. This is bla bla bla boom.
 */

Generator::Generator() :
    m_unitLength(1),
    m_nDimensions(3),
    m_idCounter(1)
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
    Vector3 offset = Vector3::createZeros();
    for(int i = 0; i < nCells; i++) {
        for(int j = 0; j < nCells; j++) {
            for(int k = 0; k < nCells; k++) {
                for(int atomi = 0; atomi < 4; atomi++) {
                    Atom* atom = new Atom(atomType);
                    Vector3 face;;
                    switch(atomi) {
                    case 0:
                        break;
                    case 1:
                        face = Vector3(sideLength / 2, sideLength / 2, 0);
                        break;
                    case 2:
                        face = Vector3(sideLength / 2, 0, sideLength / 2);
                        break;
                    case 3:
                        face = Vector3(0, sideLength / 2, sideLength / 2);
                        break;
                    }
                    Vector3 position = Vector3::createZeros();
                    position = offset + face;
                    atom->setPosition(position);
                    atom->clearDisplacement();
                    atom->setID(m_idCounter);
                    atomList.push_back(atom);
                    m_idCounter++;
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

//    cout << "Generated " << atomList.size() << " atoms in FCC structure with side length " << sideLength << endl;
//    cout << "Boundaries are " << m_lastBoundaries << endl;
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
    Vector3 totalVelocity = Vector3::createZeros();
    for(Atom* atom : atoms) {
        rowvec randVec = randn<rowvec>(m_nDimensions);
        Vector3 velocity;
        velocity[0] = randVec[0];
        velocity[1] = randVec[1];
        velocity[2] = randVec[2];
        velocity *= sqrt(temperature / atom->mass());
        totalVelocity += velocity;
        atom->setVelocity(velocity);
    }
//    cout << totalVelocity << endl;
    // Remove total linear momentum
    Vector3 velocityToRemove = totalVelocity / atoms.size();
    averageVelocity = 0;
    for(Atom* atom : atoms) {
        Vector3 newVelocity = atom->velocity() - velocityToRemove;
        atom->setVelocity(newVelocity);
        averageVelocity += dot(newVelocity, newVelocity) / atoms.size();
        totalVelocity += newVelocity;
    }
//    cout << "Boltzmann distributed velocities for " << atoms.size() << " atoms!" << endl;
//    cout << "Average velocity is " << averageVelocity << " for the temperature " << temperature << endl;
}

void Generator::uniformDistributeVelocities(double maxVelocity, vector<Atom*> atoms) {

    Vector3 totalVelocity = Vector3::createZeros();
    for(Atom* atom : atoms) {
        rowvec randVec = randu<rowvec>(m_nDimensions);
        randVec -= 0.5 * ones<rowvec>(m_nDimensions);
        Vector3 velocity(randVec);
        velocity *= 2 * maxVelocity;
        totalVelocity += velocity;
        atom->setVelocity(velocity);
    }
    // Remove total linear momentum
    Vector3 velocityToRemove = totalVelocity / atoms.size();
    totalVelocity.zeros();
    for(Atom* atom : atoms) {
        Vector3 newVelocity = atom->velocity() - velocityToRemove;
        atom->setVelocity(newVelocity);
    }
    cout << "Uniformly distributed velocities for " << atoms.size() << " atoms!" << endl;
}

const mat &Generator::lastBoundaries() const
{
    return m_lastBoundaries;
}
