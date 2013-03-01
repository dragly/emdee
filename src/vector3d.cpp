#include "vector3d.h"

Vector3::Vector3()
{
    mem_local[0] = 0;
    mem_local[1] = 0;
    mem_local[2] = 0;
}

Vector3::Vector3(double x, double y, double z)
{
    mem_local[0] = x;
    mem_local[1] = y;
    mem_local[2] = z;
}

template<class Archive>
void Vector3::serialize(Archive & ar, const unsigned int version)
{
    ar &mem_local[0];
    ar &mem_local[1];
    ar &mem_local[2];
}

std::ostream& operator<< (std::ostream &out, const Vector3 &vector)
{
    out << vector.mem_local[0] << ", " << vector.mem_local[1] << ", " << vector.mem_local[2];
    return out;
}