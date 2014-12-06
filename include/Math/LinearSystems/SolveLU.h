#ifndef SOLVE_LU_H
#define SOLVE_LU_H

/**
* Includes
**/
#include <Math/LinearSystems/LUDecomposition.h>

namespace LinearAlgebraTools {

namespace Math {

namespace LinearSystems {

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

} } }

#endif
