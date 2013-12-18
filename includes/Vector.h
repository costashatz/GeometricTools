#ifndef VECTOR_H
#define VECTOR_H

/**
* Includes
**/
#include <iostream>
#include <limits>
#include <cmath>
#include <cstring>
#include "Helper.h"
using namespace std;

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
	void Initialize(double h, Args&&... args)
	{
		if(index>=N)
			index=0;
		// add first argument to the vector data
		values[index] = h;
		index++;
		// Initialize with rest of arguments
		Initialize(forward<Args>(args)...);
	}

	/**
	* Initialize Vector from double values
	* C++11 feature - variadic templates for multiple arguments
	* Terminating Case - end of arguments
	**/
	template<typename... Args>
	void Initialize()
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
	* @params h - list of double values
	**/
	template<typename... Args>
	Vector(double h, Args&&... args)
	{
		index = 0;
		memset(values,0,N*sizeof(double));
		// use recursive initialization
		Initialize(h, args...);
	}

	/**
	* Overloading + operator
	* Addition of 2 Vectors
	* @params other - vector to perform addition with
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
	* @params other - vector to perform addition with
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
	* @params other - vector to perform subtraction with
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
	* @params other - double to multiply with
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
	* @params other - double to multiply with
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
	* @params other - double to divide with
	* @return Vector - the result
	**/
	Vector operator/(const double& other)
	{
		if(other<MinDouble)
			return (*this);
		Vector t = Vector(*this);
		for(int i=0;i<N;i++)
			t[i] /= other;
		return t;
	}

	/**
	* Overloading /= operator
	* Division with scalar (double) - if zero ignores division (returns self)
	* @params other - double to divide with
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
	double Length()
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
	double LengthSq()
	{
		double s = 0.0;
		for(int i=0;i<N;i++)
			s += values[i]*values[i];
		return s;
	}


	/**
	* Normalize Vector
	**/
	void Normalize()
	{
		double l = Length();
		if(l>std::numeric_limits<double>::epsilon())
			*this /= l;
	}

	/**
	* Overloading () operator
	* Access Vector Matlab-like
	* @params i - index to return
	* @return double - value of i-th element
	**/
	double& operator()(int i) //needs assert legal index
	{
		return values[i];
	}

	/**
	* Overloading [] operator
	* Access Vector Array-like
	* @params i - index to return
	* @return double - value of i-th element
	**/
	double& operator[](int i) //assert legal index
	{
		return values[i];
	}

	/**
	* Overloading * operator
	* Dot product with other vector
	* @params other - vector to perform dot product with
	* @return double - the dot product
	**/
	double operator*(const Vector& other)
	{
		double s=0.0;
		for(int i=0;i<N;i++)
			s += values[i]*other.values[i];
		return s;
	}

	/**
	* Print Vector Data
	**/
	void Print() // need to make it generic (overload << operator)
	{
		for(int i=0;i<N;i++)
			cout<<values[i]<<" ";
		cout<<endl;
	}

	/**
	* Destructor
	**/
	~Vector(void)
	{
	}
};


/**
* Overloading * operator
* Multiplication with scalar (double)
* Perform scalar * Vector (as opposed to Vector * scalar inside Class)
**/
template<int N>
Vector<N> operator*(const double& a, Vector<N>& b)
{
	return b*a;
}

/**
* Overloading minus (-) operator
**/
template<int N>
Vector<N> operator-(Vector<N>& b)
{
	return -1*b;
}

/**
* Defines for frequently used types
**/
#define Vector2 Vector<2>
#define Vector3 Vector<3>
#define Vector4 Vector<4>

#endif