#include <src/generator.h>

#include <src/atom.h>
#include <src/molecule.h>
/*!
 * \brief Generator::Generator has multiple functions to generate different setups of atoms.
 * Something is longer than shorter is simple. This is bla bla bla boom.
 */

Generator::Generator() :
    m_unitLength(1),
    m_boltzmannConstant(1),
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
vector<Molecule*> Generator::generateFcc(double sideLength, int nCells, AtomType atomType) {
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
                    molecule->setPosition(position);
                    molecule->clearDisplacement();
                    moleculeList.push_back(molecule);
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

    cout << "Generated " << moleculeList.size() << " molecules in FCC structure!" << endl;
    cout << "Boundaries are " << m_lastBoundaries << endl;
    return moleculeList;
}

void Generator::boltzmannDistributeVelocities(double temperature, vector<Molecule*> molecules) {

    rowvec totalVelocity = zeros<rowvec>(m_nDimensions);
    for(Molecule* molecule : molecules) {
        rowvec velocity = randn<rowvec>(m_nDimensions);
        velocity *= sqrt(m_boltzmannConstant * temperature / molecule->mass());
        totalVelocity += velocity;
        molecule->setVelocity(velocity);
    }
    // Remove total linear momentum
    rowvec velocityToRemove = totalVelocity / molecules.size();
    totalVelocity.zeros();
    double averageVelocity = 0;
    for(Molecule* molecule : molecules) {
        rowvec newVelocity = molecule->velocity() - velocityToRemove;
        molecule->setVelocity(newVelocity);
        averageVelocity += norm(newVelocity, 2) / molecules.size();
    }
    cout << "Boltzmann distributed velocities for " << molecules.size() << " molecules!" << endl;
    cout << "Average velocity is " << averageVelocity << endl;
}

void Generator::uniformDistributeVelocities(double maxVelocity, vector<Molecule*> molecules) {

    rowvec totalVelocity = zeros<rowvec>(m_nDimensions);
    for(Molecule* molecule : molecules) {
        rowvec velocity = randu<rowvec>(m_nDimensions);
        velocity -= 0.5 * ones<rowvec>(m_nDimensions);
        velocity *= 2 * maxVelocity;
        totalVelocity += velocity;
        molecule->setVelocity(velocity);
    }
    // Remove total linear momentum
    rowvec velocityToRemove = totalVelocity / molecules.size();
    totalVelocity.zeros();
    for(Molecule* molecule : molecules) {
        rowvec newVelocity = molecule->velocity() - velocityToRemove;
        molecule->setVelocity(newVelocity);
    }
    cout << "Uniformly distributed velocities for " << molecules.size() << " molecules!" << endl;
}

//void Generator::setUnitLength(double unitLength)
//{
//    m_unitLength = unitLength;
//}

const mat &Generator::lastBoundaries() const
{
    return m_lastBoundaries;
}

void Generator::setBoltzmannConstant(double boltzmannConstant) {
    m_boltzmannConstant = boltzmannConstant;
}
