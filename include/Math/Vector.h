/**
Copyright (c) 2014, Konstantinos Chatzilygeroudis
All rights reserved.
Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer
    in the documentation and/or other materials provided with the distribution.
3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived
    from this software without specific prior written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING,
BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
**/

#ifndef GEOMETRIC_TOOLS_MATH_VECTOR_H
#define GEOMETRIC_TOOLS_MATH_VECTOR_H

/**
* Includes
**/
#include <iostream>
#include <limits>
#include <cmath>
#include <cstring>
#include <Misc/Helper.h>
#include <cblas.h>

namespace GeometricTools { namespace Math {

template<unsigned int ROWS, unsigned int COLS>
class Matrix;

/**
* General Vector Class
* Supports N-dimension vectors
**/
template<unsigned int N>
class Vector
{
protected:
    // values of vector
    double values_[N];
    // helper variable for initialization
    unsigned int index_;


    /**
    * Initialize Vector from double values
    * C++11 feature - variadic templates for multiple arguments
    **/
    template<typename... Args>
    void initialize(double h, Args&&... args)
    {
        if(index_>=N)
            index_=0;
        // add first argument to the vector data
        values_[index_] = h;
        index_++;
        // Initialize with rest of arguments
        initialize(std::forward<Args>(args)...);
    }

    /**
    * Initialize Vector from double values
    * C++11 feature - variadic templates for multiple arguments
    * Terminating Case - end of arguments
    **/
    template<typename... Args>
    void initialize()
    {
        index_ = 0;
    }

public:	

    /**
    * Default Constructor
    * initialize all data to zero
    **/
    Vector()
    {
        memset(values_,0,N*sizeof(double));
    }

    /**
    * Copy Constructor
    * @param other - Vector to copy from
    **/
    Vector(const Vector& other)
    {
        memcpy(values_,other.values_,N*sizeof(double));
    }

    /**
    * Constructor
    * @param h - list of double values
    **/
    template<typename... Args>
    Vector(double h, Args&&... args)
    {
        index_ = 0;
        memset(values_,0,N*sizeof(double));
        // use recursive initialization
        initialize(h, args...);
    }

    /**
    * Overloading += operator
    * Addition of 2 Vectors
    * @param other - vector to perform addition with
    * @return Vector - self (as the result is saved there)
    **/
    Vector operator+=(const Vector& other)
    {
        cblas_daxpy(N, 1.0, &other.values_[0], 1, &values_[0], 1);
        return *this;
    }

    /**
    * Overloading -= operator
    * Subtraction of 2 Vectors
    * @param other - vector to perform subtraction with
    * @return Vector - self (as the result is saved there)
    **/
    Vector operator-=(const Vector& other)
    {
        cblas_daxpy(N, -1.0, &other.values_[0], 1, &values_[0], 1);
        return *this;
    }

    /**
    * Overloading += operator
    * Addition with scalar (double)
    * @param other - double to add with
    * @return Vector - self (result)
    **/
    Vector operator+=(const double& other)
    {
        Vector tmp;
        tmp.ones();
        cblas_daxpy(N, other, &tmp.values_[0], 1, &values_[0], 1);
        return *this;
    }

    /**
    * Overloading -= operator
    * Subtraction of scalar (double)
    * @param other - double to substract
    * @return Vector - self (result)
    **/
    Vector operator-=(const double& other)
    {
        Vector tmp;
        tmp.ones();
        cblas_daxpy(N, -other, &tmp.values_[0], 1, &values_[0], 1);
        return *this;
    }

    /**
    * Overloading *= operator
    * Multiplication with scalar (double)
    * @param other - double to multiply with
    * @return Vector - self (result)
    **/
    Vector operator*=(const double& other)
    {
        cblas_dscal(N, other, &values_[0], 1);
        return *this;
    }

    /**
    * Overloading /= operator
    * Division with scalar (double) - if zero ignores division (returns self)
    * @param other - double to divide with
    * @return Vector - self (result)
    **/
    Vector operator/=(const double& other)
    {
        if(std::abs(other) < std::numeric_limits<double>::epsilon())
            return (*this);
        cblas_dscal(N, 1.0/other, &values_[0], 1);
        return *this;
    }

    /**
    * Overloading + operator
    * Addition of 2 Vectors
    * @return Vector - the result of the addition
    **/
    Vector operator+(const Vector& lh) const
    {
        Vector t = (*this);
        t += lh;
        return t;
    }

    /**
    * Overloading - operator
    * Subtraction of 2 Vectors
    * @return Vector - the result of the subtraction
    **/
    Vector operator-(const Vector& lh) const
    {
        Vector t = (*this);
        t -= lh;
        return t;
    }

    /**
    * Overloading * operator
    * Multiplication of 2 Vectors (Dot product)
    * @return double - the result of the multiplication
    **/
    double operator*(const Vector& lh) const
    {
        double s = 0.0;
        for(unsigned int i=0;i<N;i++)
            s += (*this)[i]*lh[i];
        return s;
    }

    /**
    * Overloading + operator
    * Addition with scalar (double)
    * @param other - double to add with
    * @return Vector - the result
    **/
    Vector operator+(const double& other) const
    {
        Vector tmp = (*this);
        tmp += other;
        return tmp;
    }

    /**
    * Overloading - operator
    * Subtraction of scalar (double)
    * @param other - double to subtract
    * @return Vector - the result
    **/
    Vector operator-(const double& other) const
    {
        Vector tmp = (*this);
        tmp -= other;
        return tmp;
    }

    /**
    * Overloading * operator
    * Multiplication with scalar (double)
    * @param other - double to multiply with
    * @return Vector - the result
    **/
    Vector operator*(const double& other) const
    {
        Vector tmp = (*this);
        tmp *= other;
        return tmp;
    }

    /**
    * Overloading / operator
    * Division with scalar (double) - if zero ignores division (returns self)
    * @param other - double to divide with
    * @return Vector - the result
    **/
    Vector operator/(const double& other) const
    {
        Vector tmp = (*this);
        tmp /= other;
        return tmp;
    }

    /**
    * Overloading minus (-) operator
    **/
    Vector operator-()
    {
        return -1*(*this);
    }

    /**
    * Overloading == operator
    * @param other - Vector to compare
    * @return bool
    **/
    bool operator==(const Vector& other) const
    {
        Vector tmp = (*this)-other;
        return (tmp.lengthSq()<std::numeric_limits<double>::epsilon());
    }

    /**
    * Overloading != operator
    * @param other - Vector to compare
    * @return bool
    **/
    bool operator!=(const Vector& other) const
    {
        Vector tmp = (*this)-other;
        return (tmp.lengthSq()>=std::numeric_limits<double>::epsilon());
    }

    // !-- THIS NEEDS TO BE REMOVED --! //
    /**
    * Get pointer to values array
    * @return pointer to array
    **/
    double* data()
    {
        return values_;
    }

    /**
    * Get Length of Vector
    * @return double - length
    **/
    double length() const
    {
        return cblas_dnrm2(N, &values_[0], 1);
    }

    /**
    * Set vector to 1s
    **/
    void ones()
    {
        memset(&values_[0], 1, N*sizeof(double));
    }


    /**
    * Get LengthSq of Vector
    * @return double - length squared
    **/
    double lengthSq() const
    {
        double s = cblas_dnrm2(N, &values_[0], 1);
        return s*s;
    }


    /**
    * Normalize Vector
    **/
    void normalize()
    {
        double l = length();
        if(l>std::numeric_limits<double>::epsilon())
            *this /= l;
    }

    /**
    * Normalized Vector
    **/
    Vector normalized()
    {
        Vector tmp = Vector(*this);
        tmp.normalize();
        return tmp;
    }

    /**
    * Overloading () operator
    * Access Vector Matlab-like
    * @param i - index to return
    * @return double - value of i-th element
    **/
    double& operator()(int i) //needs assert legal index
    {
        return values_[i];
    }

    /**
    * Overloading () operator - const version
    * Access Vector Matlab-like
    * @param i - index to return
    * @return double - value of i-th element
    **/
    double operator()(int i) const //needs assert legal index
    {
        return values_[i];
    }

    /**
    * Overloading [] operator
    * Access Vector Array-like
    * @param i - index to return
    * @return double - value of i-th element
    **/
    double& operator[](int i) //assert legal index
    {
        return values_[i];
    }

    /**
    * Overloading [] operator - const version
    * Access Vector Array-like
    * @param i - index to return
    * @return double - value of i-th element
    **/
    double operator[](int i) const //assert legal index
    {
        return values_[i];
    }

    /**
    * Get ith unit vector with size N
    * @param i
    * @return Vector - unit vector
    **/
    static Vector e(const unsigned int& i)
    {
        Vector k;
        k[i] = 1;
        return k;
    }

    /**
    * Destructor
    **/
    ~Vector(void)
    {
    }

    template<unsigned int ROWS, unsigned int COLS>
    friend class Matrix;
};

/**
* Cross Product
* Applies only to 3D Vectors
* @return Vector - the result of the cross product
**/
inline Vector<3> cross(const Vector<3>& rh, const Vector<3>& lh)
{
    return Vector<3>(rh[1]*lh[2]-rh[2]*lh[1], rh[2]*lh[0]-rh[0]*lh[2], rh[0]*lh[1]-rh[1]*lh[0]);
}

/**
* Projection of a to e
* @return Vector - the result of the cross product
**/
template <unsigned int N>
Vector<N> projection(const Vector<N>& a, const Vector<N>& e)
{
    return ((e*a)*e/(e*e));
}

/**
* Overloading << operator
* "print" vector to stream
* @param ostream - stream to print the vector to
* @param Vector - vector to print
**/
template <unsigned int N>
std::ostream& operator<<(std::ostream& os, const Vector<N>& obj)
{
    for(unsigned int i=0;i<N-1;i++)
        os<<obj[i]<<" ";
    os<<obj[N-1];
    return os;
}

/**
* Overloading >> operator
* read vector from stream
* @param istream - stream to read the vector from
* @param Vector - vector to read
**/
template <unsigned int N>
std::istream& operator>>(std::istream& in, Vector<N>& obj)
{
    for(unsigned int i=0;i<N;i++)
        in>>obj[i];
    return in;
}

/**
* Overloading / operator
* Division with scalar (double)
* Perform Vector / scalar (as opposed to Vector / scalar inside Class)
**/
template <unsigned int N>
Vector<N> operator/(const double& a, const Vector<N>& b)
{
    if(std::abs(a) < std::numeric_limits<double>::epsilon())
        return Vector<N>();
    Vector<N> tmp = Vector<N>(b);
    for(unsigned int i=0;i<N;i++) {
        if(std::abs(tmp[i]) < std::numeric_limits<double>::epsilon())
            return b;
        tmp[i] = a/tmp[i];
    }
    return tmp;
}

/**
* Overloading * operator
* Multiplication with scalar (double)
* Perform scalar * Vector (as opposed to Vector * scalar inside Class)
**/
template <unsigned int N>
Vector<N> operator*(const double& a, const Vector<N>& b)
{
    return b*a;
}

/**
* Typedefs for frequently used types
**/
typedef Vector<2> Vector2;
typedef Vector<3> Vector3;
typedef Vector<4> Vector4;

} }

#endif
