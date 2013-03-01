#ifndef VECTOR3D_H
#define VECTOR3D_H

#include <iostream>
#include <boost/serialization/serialization.hpp>

class Vector3
{
public:
    Vector3();
    Vector3(double x, double y, double z);

    double x() const;
    double y() const;
    double z() const;

    double& operator()(const int component);
    double& operator[](const int component);
    friend std::ostream& operator<< (std::ostream &out, const Vector3 &vector);
    friend Vector3 operator+ (const Vector3 &vector1, const Vector3 &vector2);
    friend Vector3 operator- (const Vector3 &vector1, const Vector3 &vector2);
    friend double operator* (const Vector3 &vector1, const Vector3 &vector2);
    friend bool operator== (const Vector3 &vector1, const Vector3 &vector2);
    friend bool operator!= (const Vector3 &vector1, const Vector3 &vector2);
protected:
    double mem_local[3];

private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
};


inline double& Vector3::operator[](const int component)
{
    return mem_local[component];
}

inline double& Vector3::operator()(const int component)
{
    return mem_local[component];
}

inline Vector3 operator+ (const Vector3 &vector1, const Vector3 &vector2)
{
    return Vector3(vector1.mem_local[0] + vector2.mem_local[0], vector1.mem_local[1] + vector2.mem_local[1], vector1.mem_local[2] + vector2.mem_local[2]);
}

inline Vector3 operator- (const Vector3 &vector1, const Vector3 &vector2)
{
    return Vector3(vector1.mem_local[0] - vector2.mem_local[0], vector1.mem_local[1] - vector2.mem_local[1], vector1.mem_local[2] - vector2.mem_local[2]);
}

/*!
 * \brief Vector3::operator * returns the dot product of the two vectors.
 * \param vector1
 * \param vector2
 * \return
 */
inline double operator* (const Vector3 &vector1, const Vector3 &vector2)
{
    return (vector1.mem_local[0] * vector2.mem_local[0] + vector1.mem_local[1] * vector2.mem_local[1] + vector1.mem_local[2] * vector2.mem_local[2]);
}


inline bool operator== (const Vector3 &vector1, const Vector3 &vector2)
{
    return ((vector1.mem_local[0] == vector2.mem_local[0]) &&
            (vector1.mem_local[1] == vector2.mem_local[1]) &&
            (vector1.mem_local[2] == vector2.mem_local[2]));
}

inline bool operator != (const Vector3 &vector1, const Vector3 &vector2)
{
    return !(vector1 == vector2);
}

double Vector3::x() const {
    return mem_local[0];
}
double Vector3::y() const {
    return mem_local[1];
}
double Vector3::z() const {
    return mem_local[2];
}

#endif // VECTOR3D_H
