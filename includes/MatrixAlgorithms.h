#ifndef MATRIX_ALGO_H
#define MATRIX_ALGO_H

#include "Matrix.h"

/**
* Matrix (Linear Algebra) Algorithms
* So far:
*	*Solve Linear System using Gauss Elimination
*	*LU Decomposition
*	*Solve Linear System using LU Decomposition
**/


/**
* Solve Linear System using Gauss Elimination
* @param A - parameter Matrix
* @param B - constant Vector
* @return Vector<D> - solution
**/
template<unsigned int D>
Vector<D> solveGauss(const Matrix<D,D>& A, const Vector<D>& B)
{
	Matrix<D,D> a = A;
	Vector<D> b = B;
	Vector<D> res;
	for(unsigned int i=0;i<D-1;i++)
	{
		int m = 0;
		double temp1 = 0.0;
		for(unsigned int j=i;j<D;j++)
		{
			double temp2 = 0.0;
			for(unsigned int k=i+1;k<D;k++)
				if(abs(a(j,k))>temp2)
					temp2 = abs(a(j,k));
			if(abs(a(j,i))>(temp1*temp2))
			{
				temp1 = abs(a(j,i))/temp2;
				m = j;
			}
		}
		if(m!=i)
		{
			a.SwapRows(i,m);
			double temp = b(i);
			b(i) = b(m);
			b(m) = temp;
		}
		
		if(abs(a(i,i))<std::numeric_limits<double>::epsilon())
			return res;

		for(unsigned int k=i+1;k<D;k++)
		{
			double r = a(k,i)/a(i,i);
			for(unsigned int c=0;c<D;c++)
			{
				a(k,c) -= r*a(i,c);
			}
			b(k) -= r*b(i);
		}
	}
	for(int k=D-1;k>=0;k--)
	{
		double S = 0.0;
		for(unsigned int j=k+1;j<D;j++)
		{
			S += a(k,j)*res(j);
		}
		res(k) = 1.0/a(k,k)*(b(k)-S);
	}
	return res;
}


/**
* Perform LU Decomposition
* @param a - Matrix to decompose
* @param L - L Matrix --
* @param U - U Matrix  |
* @param P - P Matrix  |
*					   --> Pass by reference (returns)
**/
template<unsigned int D>
void LUDecomposition(const Matrix<D,D>& a, Matrix<D,D>& L, Matrix<D,D>& U, Matrix<D,D>& P)
{
	L.Identity();
	P.Identity();
	U = a;
	for(int k=0;k<D-1;k++)
	{
		double m = U(k,k);
		int newJ = k;
		for(int i=k;i<D;i++)
		{
			if(abs(U(i,k))>m)
			{
				m = abs(U(i,k));
				newJ = i;
			}
		}
		if(newJ!=k)
		{
			for(int j=k;j<D;j++)
			{
				double t = U(k,j);
				U(k,j) = U(newJ,j);
				U(newJ,j) = t;
			}
			for(int j=0;j<=k-1;j++)
			{
				double t = L(k,j);
				L(k,j) = L(newJ,j);
				L(newJ,j) = t;
			}
			P.SwapRows(k,newJ);
		}

		for(int j=k+1;j<D;j++)
		{
			L(j,k) = U(j,k)/U(k,k);
			for(int i=k;i<D;i++)
			{
				U(j,i) -= L(j,k)*U(k,i);
			}
		}
	}
}

/**
* Solve Linear System using LU Decomposition
* @param A - parameter Matrix
* @param B - constant Vector
* @return Vector<D> - solution
**/
template<unsigned int D>
Vector<D> solveLU(const Matrix<D,D>& A, const Vector<D>& B)
{
	Matrix<D,D> L,U,P;
	LUDecomposition(A,L,U,P);

	Vector<D> res, y;
	Vector<D> b = P*B;


	for(int k=0;k<D;k++)
	{
		double S = 0.0;
		for(int j=k-1;j>=0;j--)
		{
			S += L(k,j)*y(j);
		}
		y(k) = 1.0/L(k,k)*(b(k)-S);
	}

	for(int k=D-1;k>=0;k--)
	{
		double S = 0.0;
		for(unsigned int j=k+1;j<D;j++)
		{
			S += U(k,j)*res(j);
		}
		res(k) = 1.0/U(k,k)*(y(k)-S);
	}

	return res;
}

#endif