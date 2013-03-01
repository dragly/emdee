#ifndef VECTOR3D_H
#define VECTOR3D_H

#include <iostream>
#include <boost/serialization/serialization.hpp>
#include <armadillo>

class Vector3
{
public:
    Vector3();
    Vector3(double x, double y, double z);
    Vector3(arma::rowvec armaVector);

    // Operations
    void zeros();

    // Accessors
    double x() const;
    double y() const;
    double z() const;

    uint size() const;

    // Operator overloading
    double& operator()(const int component);
    double& operator[](const int component);
    double operator()(const int component) const;
    double operator[](const int component) const;
    friend std::ostream& operator<< (std::ostream &out, const Vector3 &vector);
    friend Vector3 operator+ (const Vector3 &vector1, const Vector3 &vector2);
    friend Vector3 operator- (const Vector3 &vector1, const Vector3 &vector2);
    friend Vector3 operator+(const Vector3 &vector1);
    friend Vector3 operator-(const Vector3 &vector1);
    Vector3& operator= (const Vector3 &vector2);
    Vector3& operator+= (const Vector3 &vector2);
    Vector3& operator-= (const Vector3 &vector2);
    friend double operator* (const Vector3 &vector1, const Vector3 &vector2);
    friend bool operator== (const Vector3 &vector1, const Vector3 &vector2);
    friend bool operator!= (const Vector3 &vector1, const Vector3 &vector2);

    // Operators with doubles
    friend Vector3 operator*(const Vector3 &vector1, double value);
    friend Vector3 operator*(double value, const Vector3 &vector1);
    friend Vector3 operator/(const Vector3 &vector1, double value);
    Vector3& operator*= (double value);
    Vector3& operator/= (double value);

    // Statics
    static Vector3 ones();
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


inline double Vector3::operator[](const int component) const
{
    return mem_local[component];
}

inline double Vector3::operator()(const int component) const
{
    return mem_local[component];
}

inline void Vector3::zeros() {
    mem_local[0] = 0;
    mem_local[1] = 0;
    mem_local[2] = 0;
}

inline uint Vector3::size() const {
    return 3;
}


inline Vector3 operator+ (const Vector3 &vector1, const Vector3 &vector2)
{
    return Vector3(vector1.mem_local[0] + vector2.mem_local[0], vector1.mem_local[1] + vector2.mem_local[1], vector1.mem_local[2] + vector2.mem_local[2]);
}

inline Vector3 operator- (const Vector3 &vector1, const Vector3 &vector2)
{
    return Vector3(vector1.mem_local[0] - vector2.mem_local[0], vector1.mem_local[1] - vector2.mem_local[1], vector1.mem_local[2] - vector2.mem_local[2]);
}

inline Vector3& Vector3::operator= (const Vector3 &vector2)
{
    this->mem_local[0] = vector2.mem_local[0];
    this->mem_local[1] = vector2.mem_local[1];
    this->mem_local[2] = vector2.mem_local[2];
    return *this;
}

inline Vector3& Vector3::operator+= (const Vector3 &vector2)
{
    this->mem_local[0] += vector2.mem_local[0];
    this->mem_local[1] += vector2.mem_local[1];
    this->mem_local[2] += vector2.mem_local[2];
    return *this;
}

inline Vector3& Vector3::operator*= (double value)
{
    this->mem_local[0] *= value;
    this->mem_local[1] *= value;
    this->mem_local[2] *= value;
    return *this;
}

inline Vector3& Vector3::operator/= (double value)
{
    this->mem_local[0] /= value;
    this->mem_local[1] /= value;
    this->mem_local[2] /= value;
    return *this;
}

inline Vector3 operator*(const Vector3 &vector1, double value)
{
    return Vector3(vector1.mem_local[0]*value, vector1.mem_local[1]*value, vector1.mem_local[2]*value);
}

inline Vector3 operator*(double value, const Vector3 &vector1)
{
    return Vector3(vector1.mem_local[0]*value, vector1.mem_local[1]*value, vector1.mem_local[2]*value);
}

inline Vector3 operator/(const Vector3 &vector1, double value)
{
    return Vector3(vector1.mem_local[0]/value, vector1.mem_local[1]/value, vector1.mem_local[2]/value);
}

inline Vector3& Vector3::operator-= (const Vector3 &vector2)
{
    this->mem_local[0] -= vector2.mem_local[0];
    this->mem_local[1] -= vector2.mem_local[1];
    this->mem_local[2] -= vector2.mem_local[2];
    return *this;
}

inline Vector3 operator+(const Vector3 &vector1)
{
    return Vector3(vector1.mem_local[0], vector1.mem_local[1], vector1.mem_local[2]);
}

inline Vector3 operator-(const Vector3 &vector1)
{
    return Vector3(-vector1.mem_local[0], -vector1.mem_local[1], -vector1.mem_local[2]);
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

inline double Vector3::x() const {
    return mem_local[0];
}
inline double Vector3::y() const {
    return mem_local[1];
}
inline double Vector3::z() const {
    return mem_local[2];
}

#endif // VECTOR3D_H
