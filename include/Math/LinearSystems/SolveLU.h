#ifndef SOLVE_LU_H
#define SOLVE_LU_H

/**
* Includes
**/
#include <Math/Matrix.h>
#include <lapacke.h>

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
    Matrix<D,D> U = A;
    int* ipiv = new int[D];
    LAPACKE_dgetrf(LAPACK_ROW_MAJOR, D, D, U.values, D, ipiv);

    Vector<D> res = B;
    LAPACKE_dgetrs(LAPACK_ROW_MAJOR, 'N', D, 1, U.values, D, ipiv, res.values, 1);

    return res;
}

} } }

#endif
