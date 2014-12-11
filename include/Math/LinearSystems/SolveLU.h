#ifndef SOLVE_LU_H
#define SOLVE_LU_H

/**
* Includes
**/
#include <Math/Matrix.h>

extern "C" int dgetrf_(int *m, int *n, double *a, int * lda, int *ipiv, int *info);
extern "C" int dgetrs_(char *trans, int *n, int *nrhs,
                       double *a, int *lda, int *ipiv, double *b, int * ldb, int *info);

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
    int dim = int(D), info, nhrs = 1;
    char trans = 'T';
    Matrix<D,D> U = A;
    int* ipiv = new int[dim];

    dgetrf_(&dim, &dim, U.data(), &dim, ipiv, &info);

    Vector<D> res = B;

    dgetrs_(&trans, &dim, &nhrs, U.data(), &dim, ipiv, res.data(), &dim, &info);

    return res;
}

} } }

#endif
