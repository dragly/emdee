#include "fannforce.h"

#include <doublefann.h>
#include <atom.h>

FannForce::FannForce() :
    m_ann(0),
    m_hasWarnedAboutMissingNetwork(false)
{
}

void FannForce::loadNetwork(std::string fileName)
{
    m_ann = fann_create_from_file(fileName.c_str());
}

//void FannForce::calculateAndApplyForce(Atom *atom1, Atom *atom2, Atom *atom3)
//{
//    (void)atom1;
//    (void)atom2;
//    (void)atom3;
//    if(!m_ann) {
//        if(!m_hasWarnedAboutMissingNetwork) {
//            cerr << "FANN network not loaded. Cannot apply force." << endl;
//            m_hasWarnedAboutMissingNetwork = true;
//        }
//        return;
//    }
//}

//void FannForce::calculateAndApplyForce(Atom *atom1, Atom *atom2, Atom *atom3, const Vector3 &atom2Offset, const Vector3 &atom3Offset)
//{
//    (void)atom1;
//    (void)atom2;
//    (void)atom3;
//    (void)atom2Offset;
//    (void)atom3Offset;
//    if(!m_ann) {
//        if(!m_hasWarnedAboutMissingNetwork) {
//            cerr << "FANN network not loaded. Cannot apply force." << endl;
//            m_hasWarnedAboutMissingNetwork = true;
//        }
//        return;
//    }
//}

void FannForce::calculateAndApplyForce(Atom *atom1, Atom *atom2, Atom *atom3)
{
    atom1->addForce(0, 0.1);
    atom2->addForce(0, 0.1);
    atom3->addForce(0, 0.1);
}

void FannForce::calculateAndApplyForce(Atom *atom1, Atom *atom2, Atom *atom3, const Vector3 &atom2Offset, const Vector3 &atom3Offset)
{

}
