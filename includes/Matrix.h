#ifndef MATRIX_H
#define MATRIX_H

/**
* Includes and Defines
**/
#include <iostream>
#include <cassert>
#include "Vector.h"
#include "CStrassen.h"
using namespace std;

#define MAX(a,b) a>b?a:b

/**
* General Matrix Class
* Supports NxM-dimension matrices
**/
template<unsigned int ROWS, unsigned int COLS>
class Matrix
{
protected:
	// matrix data
	double values[COLS*ROWS+1];
	// helper variable for initialization
	unsigned int index;

	/**
	* Initialize Matrix from double values
	* C++11 feature - variadic templates for multiple arguments
	**/
	template<typename... Args>
	void Initialize(double h, Args&&... args)
	{
		if(index>=(COLS*ROWS))
			index=0;
		unsigned int i = index/COLS;
		unsigned int j = index%COLS;
		// add first argument to the matrix data
		values[i*COLS+j] = h;
		index++;
		// Initialize with rest of arguments
		Initialize(forward<Args>(args)...);
	}

	/**
	* Initialize Matrix from double values
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
	Matrix()
	{
		memset(values,0,COLS*ROWS*sizeof(double));
	}

	/**
	* Copy Constructor
	* @param other - Matrix to copy from
	**/
	Matrix(const Matrix& other)
	{
		memcpy(values,other.values,COLS*ROWS*sizeof(double));
	}

	/**
	* Constructor
	* @params h - list of double values
	**/
	template<typename... Args>
	Matrix(double h, Args&&... args)
	{
		index = 0;
		memset(values,0,COLS*ROWS*sizeof(double));
		Initialize(h, args...);
	}

	/**
	* Make Identity
	* Applies only to square matrices
	* If not square, nothing happens
	**/
	void Identity()
	{
		if(COLS!=ROWS)
			return;
		memset(values,0,COLS*ROWS*sizeof(double));
		for(unsigned int i=0;i<COLS;i++)
			(*this)(i,i) = 1.0;
	}

	/**
	* Make Zero
	**/
	void Zero()
	{
		memset(values,0,COLS*ROWS*sizeof(double));
	}

	/**
	* Overloading () operator
	* Access Matrix Matlab-like
	* @params i - row index to return
	* @params j - column index to return
	* @return double - value of i-th element
	**/
	double& operator()(unsigned int i, unsigned int j) //assert legal index
	{
		return values[i*COLS+j];
	}

	/**
	* Overloading + operator
	* Addition of 2 Matrices
	* @params other - Matrix to perform addition with
	* @return Matrix - the result of the addition
	**/
	Matrix operator+(const Matrix& other)
	{
		unsigned int N = COLS*ROWS;
		Matrix t = Matrix(*this);
		for(unsigned int i=0;i<N;i++)
			t.values[i] += other.values[i];
		return t;
	}

	Matrix operator+=(const Matrix& other)
	{
		*this = *this+other;
		return *this;
	}

	/**
	* Overloading - operator
	* Addition of 2 Matrices
	* @params other - Matrix to perform subtraction with
	* @return Matrix - the result
	**/
	Matrix operator-(const Matrix& other)
	{
		unsigned int N = COLS*ROWS;
		Matrix t = Matrix(*this);
		for(unsigned int i=0;i<N;i++)
			t.values[i] -= other.values[i];
		return t;
	}

	Matrix operator-=(const Matrix& other)
	{
		*this = *this-other;
		return *this;
	}

	/**
	* Overloading * operator
	* Multiplication with scalar (double)
	* @params other - double to multiply with
	* @return Matrix - the result
	**/
	Matrix operator*(const double& other)
	{
		unsigned int N = COLS*ROWS;
		Matrix t = Matrix(*this);
		for(unsigned int i=0;i<N;i++)
			t.values[i] *= other;
		return t;
	}

	Matrix operator*=(const double& other)
	{
		*this = *this*other;
		return *this;
	}

	/**
	* Overloading / operator
	* Division with scalar (double) - if zero ignores division (returns self)
	* @params other - double to divide with
	* @return Matrix - the result
	**/
	Matrix operator/(const double& other)
	{
		if(other<MinDouble)
			return (*this);
		unsigned int N = COLS*ROWS;
		Matrix t = Matrix(*this);
		for(unsigned int i=0;i<N;i++)
			t.values[i] /= other;
		return t;
	}

	Matrix operator/=(const double& other)
	{
		*this = *this/other;
		return *this;
	}

	/**
	* Get Transpose of Matrix
	* @return Matrix - the transposed matrix
	**/
	Matrix<COLS,ROWS> Transpose()
	{
		Matrix<COLS,ROWS> r;
		for(unsigned int i=0;i<ROWS;i++)
		{
			for(unsigned int j=0;j<COLS;j++)
				r(j,i) = (*this)(i,j);
		}
		return r;
	}

	Matrix<COLS,ROWS> operator~()
	{
		Matrix<COLS,ROWS> r = this->Transpose();
		return r;
	}

	/**
	* Print Matrix Data
	**/
	void Print() // need to make it generic (overload << operator)
	{
		for(unsigned int i=0;i<ROWS;i++)
		{
			for(unsigned int j=0;j<COLS;j++)
				cout<<(*this)(i,j)<<" ";
			cout<<endl;
		}
	}

	/**
	* Get Norm of Matrix
	* @return double - the L2 Norm of the Matrix
	**/
	double Norm()
	{
		double s = 0.;
		for(unsigned int i=0;i<ROWS;i++)
		{
			for(unsigned int j=0;j<COLS;j++)
				s += (*this)(i,j)*(*this)(i,j);
		}
		return sqrt(s);
	}

	/**
	* Swap Row i with Row j
	* @param i - row index i
	* @param j - row index j
	**/
	void SwapRows(unsigned int i, unsigned int j)
	{
		Vector<COLS> t;
		for(unsigned k=0;k<COLS;k++)
			t(k) = (*this)(i,k);
		for(unsigned k=0;k<COLS;k++)
		{
			(*this)(i,k) = (*this)(j,k);
			(*this)(j,k) = t(k);
		}
	}

	/**
	* Swap Column i with Column j
	* @param i - column index i
	* @param j - column index j
	**/
	void SwapCols(unsigned int i, unsigned int j)
	{
		Vector<ROWS> t;
		for(unsigned k=0;k<ROWS;k++)
			t(k) = (*this)(k,i);
		for(unsigned k=0;k<COLS;k++)
		{
			(*this)(k,i) = (*this)(k,j);
			(*this)(k,j) = t(k);
		}
	}

	/**
	* Get # rows of Matrix
	* @return unsigned int - # rows
	**/
	unsigned int rows()
	{
		return ROWS;
	}

	/**
	* Get # columns of Matrix
	* @return unsigned int - # columns
	**/
	unsigned int cols()
	{
		return COLS;
	}
};

/**
* Overloading * operator
* Multiplication with scalar (double)
* Perform scalar * Matrix (as opposed to Matrix * scalar inside Class)
**/
template<unsigned int N, unsigned int M>
Matrix<N,M> operator*(const double& a, Matrix<N,M>& b)
{
	return b*a;
}

/**
* Overloading minus (-) operator
**/
template<unsigned int N, unsigned int M>
Matrix<N,M> operator-(Matrix<N,M>& b)
{
	return -1*b;
}

/**
* Overloading * operator
* Multiplication with Vector
**/
template<unsigned int C1, unsigned int K>
Vector<C1> operator*(const Matrix<C1,K>& r1, const Vector<K>& r2)
{
	Matrix<K,1> b;
	Vector<K> v = r2;
	for(unsigned int i=0;i<K;i++)
		b(i,0) = v[i];
	Matrix<C1,1> c;
	c = r1*b;
	Vector<C1> res;
	for(unsigned int i=0;i<K;i++)
		res[i] = c(i,0);
	return res;
}

/**
* Get Minor of the Matrix
**/
template<unsigned int ROWS, unsigned int COLS>
Matrix<ROWS-1,COLS-1> Minor(Matrix<ROWS,COLS> m, unsigned int i, unsigned int j)
{
	Matrix<ROWS-1,COLS-1> res;
	unsigned int rk=0, ck=0;
	for(unsigned int r=0;r<ROWS;r++)
	{
		if(r==i)
			continue;
		ck = 0;
		for(unsigned int c=0;c<COLS;c++)
		{
			if(c==j)
				continue;
			res(rk,ck) = m(r,c);
			ck++;
		}
		rk++;
	}
	return res;
}

/**
* Get Determinant of the Matrix
**/
template<unsigned int ROWS, unsigned int COLS>
double Determinant(Matrix<ROWS,COLS> m)
{
	if(COLS!=ROWS || COLS<2)
		return 0.0;
	if(COLS==2)
	{
		double x = (m(0,0)*m(1,1)-m(1,0)*m(0,1));
		return x;
	}
	double s = 0.0;
	const int R = ROWS-1;
	const int C = COLS-1;
	Matrix<R,C> temp;
	for(unsigned int i=0;i<COLS;i++)
	{
		temp = Minor(m,0,i);
		s += pow(-1,i)*m(0,i)*Determinant(temp);
	}
	return s;
}

/**
* Get Determinant of the Matrix
* Base Case - 1x1 Matrix - Scalar value
**/
template<>
double Determinant(Matrix<1,1> m)
{
	return 0.0;
}

/**
* Overloading * operator
* Multiplication with Matrix
**/
template<unsigned int C1, unsigned int K, unsigned int R2>
Matrix<C1,R2> operator*(const Matrix<C1,K>& r1, const Matrix<K,R2>& r2)
{
	Matrix<C1,K> r = r1;
	Matrix<K,R2> rb = r2;
	double** a, **b, **c;
	double D1=(ceil(log(K)/log(2.))), D2=(ceil(log(C1)/log(2.))),D3=(ceil(log(R2)/log(2.)));
	unsigned int DIM = MAX(D1,MAX(D2,D3));
	DIM = (unsigned int)pow(2,DIM);
	a = allocate_real_matrix(DIM,0);
	b = allocate_real_matrix(DIM,0);
	c = allocate_real_matrix(DIM,0);
	for(unsigned int i=0;i<C1;i++)
	{
		for(unsigned int j=0;j<K;j++)
		{
			a[i][j] = r(i,j);
		}
	}
	for(unsigned int i=0;i<K;i++)
	{
		for(unsigned int j=0;j<R2;j++)
		{
			b[i][j] = rb(i,j);
		}
	}
	strassen(a,b,c,DIM);
	a = free_real_matrix(a,DIM);
	b = free_real_matrix(b,DIM);
	Matrix<C1,R2> res;
	for(unsigned int i=0;i<C1;i++)
	{
		for(unsigned int j=0;j<R2;j++)
		{
			res(i,j) = c[i][j];
		}
	}
	c = free_real_matrix(c,DIM);
	return res;
}

#endif