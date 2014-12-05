#ifndef VECTOR_H
#define VECTOR_H

/**
* Includes
**/
#include <iostream>
#include <limits>
#include <cmath>
#include <cstring>
#include <Misc/Helper.h>

namespace LinearAlgebraTools { namespace Math {

/**
* General Vector Class
* Supports N-dimension vectors
**/
template<unsigned int N=2>
class Vector
{
protected:
    // values of vector
    double values[N];
    // helper variable for initialization
    unsigned int index;


    /**
    * Initialize Vector from double values
    * C++11 feature - variadic templates for multiple arguments
    **/
    template<typename... Args>
    void initialize(double h, Args&&... args)
    {
        if(index>=N)
            index=0;
        // add first argument to the vector data
        values[index] = h;
        index++;
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
        index = 0;
    }

public:	

    /**
    * Default Constructor
    * initialize all data to zero
    **/
    Vector()
    {
        memset(values,0,N*sizeof(double));
    }

    /**
    * Copy Constructor
    * @param other - Vector to copy from
    **/
    Vector(const Vector& other)
    {
        memcpy(values,other.values,N*sizeof(double));
    }

    /**
    * Constructor
    * @param h - list of double values
    **/
    template<typename... Args>
    Vector(double h, Args&&... args)
    {
        index = 0;
        memset(values,0,N*sizeof(double));
        // use recursive initialization
        initialize(h, args...);
    }

    /**
    * Overloading + operator
    * Addition of 2 Vectors
    * @param other - vector to perform addition with
    * @return Vector - the result of the addition
    **/
    Vector operator+(const Vector& other)
    {
        Vector t = Vector(*this);
        for(int i=0;i<N;i++)
            t[i] += other.values[i];
        return t;
    }

    /**
    * Overloading += operator
    * Addition of 2 Vectors
    * @param other - vector to perform addition with
    * @return Vector - self (as the result is saved there)
    **/
    Vector operator+=(const Vector& other)
    {
        *this = *this+other;
        return *this;
    }

    /**
    * Overloading - operator
    * Subtraction of 2 Vectors
    * @params other - vector to perform subtraction with
    * @return Vector - the result of the subtraction
    **/
    Vector operator-(const Vector& other)
    {
        Vector t = Vector(*this);
        for(int i=0;i<N;i++)
            t[i] -= other.values[i];
        return t;
    }

    /**
    * Overloading -= operator
    * Subtraction of 2 Vectors
    * @param other - vector to perform subtraction with
    * @return Vector - self (as the result is saved there)
    **/
    Vector operator-=(const Vector& other)
    {
        *this = *this-other;
        return *this;
    }

    /**
    * Overloading * operator
    * Multiplication with scalar (double)
    * @param other - double to multiply with
    * @return Vector - the result
    **/
    Vector operator*(const double& other)
    {
        Vector t = Vector(*this);
        for(int i=0;i<N;i++)
            t[i] *= other;
        return t;
    }

    /**
    * Overloading *= operator
    * Multiplication with scalar (double)
    * @param other - double to multiply with
    * @return Vector - self (result)
    **/
    Vector operator*=(const double& other)
    {
        *this = *this*other;
        return *this;
    }

    /**
    * Overloading / operator
    * Division with scalar (double) - if zero ignores division (returns self)
    * @param other - double to divide with
    * @return Vector - the result
    **/
    Vector operator/(const double& other)
    {
        if(std::abs(other) < std::numeric_limits<double>::epsilon())
            return (*this);
        Vector t = Vector(*this);
        for(int i=0;i<N;i++)
            t[i] /= other;
        return t;
    }

    /**
    * Overloading /= operator
    * Division with scalar (double) - if zero ignores division (returns self)
    * @param other - double to divide with
    * @return Vector - self (result)
    **/
    Vector operator/=(const double& other)
    {
        *this = *this/other;
        return *this;
    }

    /**
    * Get Length of Vector
    * @return double - length
    **/
    double length()
    {
        double s = 0.0;
        for(int i=0;i<N;i++)
            s += values[i]*values[i];
        return sqrt(s);
    }


    /**
    * Get LengthSq of Vector
    * @return double - length squared
    **/
    double lengthSq()
    {
        double s = 0.0;
        for(int i=0;i<N;i++)
            s += values[i]*values[i];
        return s;
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
        return values[i];
    }

    /**
    * Overloading [] operator
    * Access Vector Array-like
    * @param i - index to return
    * @return double - value of i-th element
    **/
    double& operator[](int i) //assert legal index
    {
        return values[i];
    }

    /**
    * Overloading * operator
    * Dot product with other vector
    * @param other - vector to perform dot product with
    * @return double - the dot product
    **/
    double operator*(const Vector& other)
    {
        double s=0.0;
        for(int i=0;i<N;i++)
            s += values[i]*other.values[i];
        return s;
    }

    template <unsigned int U>
    friend std::ostream& operator<<(std::ostream& os, const Vector<U>& obj);
    template <unsigned int U>
    friend std::istream& operator>>(std::istream& in, Vector<U>& obj);

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
};

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
* Overloading minus (-) operator
**/
template <unsigned int N>
Vector<N> operator-(const Vector<N>& b)
{
    return -1*b;
}

/**
* Typedefs for frequently used types
**/
typedef Vector<2> Vector2;
typedef Vector<3> Vector3;
typedef Vector<4> Vector4;

} }

#endif
