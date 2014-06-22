#ifdef USE_MPI
#include <boost/mpi.hpp>
#endif

#include "vector3.h"

/*!
 * \class Vector3
 * \brief The Vector3 class is a simple and performance-centric 3D vector class.
 */

Vector3 Vector3::m_zeroVectorPointer = Vector3(0,0,0);

std::ostream& operator<< (std::ostream &out, const Vector3 &vector)
{
    out << vector.mem_local[0] << ", " << vector.mem_local[1] << ", " << vector.mem_local[2];
    return out;
}

Vector3 Vector3::ones()
{
    return Vector3(1,1,1);
}


Vector3 Vector3::createZeros()
{
    return Vector3(0,0,0);
}
