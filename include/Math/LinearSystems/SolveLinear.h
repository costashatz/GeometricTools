#ifndef SOLVE_LINEAR_H
#define SOLVE_LINEAR_H

/**
* Includes
**/
#include <Math/Matrix.h>

extern "C" void dgesv_( int* n, int* nrhs, double* a, int* lda, int* ipiv,
                double* b, int* ldb, int* info );

namespace LinearAlgebraTools {

namespace Math {

namespace LinearSystems {

/**
* Solve Linear System
* @param A - parameter Matrix
* @param B - constant Vector
* @return Vector<D> - solution
**/
template<unsigned int D>
Vector<D> solveLinear(const Matrix<D,D>& A, const Vector<D>& B)
{
    int dim = int(D), nrhs = 1, info;
    Matrix<D,D> U = A;
    int* ipiv = new int[D];

    Vector<D> res = B;
    dgesv_( &dim, &nrhs, U.transpose().data(), &dim, ipiv, res.data(), &dim, &info );


    return res;
}

} } }

#endif
