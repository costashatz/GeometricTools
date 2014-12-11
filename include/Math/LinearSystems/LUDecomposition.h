#ifndef LU_DECOMPOSITION
#define LU_DECOMPOSITION

/**
* Includes
**/
#include <Math/Matrix.h>

extern "C" int dgetrf_(int *m, int *n, double *a, int * lda, int *ipiv, int *info);

namespace LinearAlgebraTools {

namespace Math {

template<unsigned int ROWS, unsigned int COLS>
class Matrix;

template<unsigned int N>
class Vector;

namespace LinearSystems {

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
    L.identity();
    P.identity();
    U = a;
    int* ipiv = new int[D];
    int dim = int(D), info;
    dgetrf_(&dim, &dim, U.data(), &dim, ipiv, &info);
    for(unsigned int i=0;i<D;i++) {
        P.swapCols(i, ipiv[i]-1);
        for(unsigned int j=0;j<i;j++) {
            L(i,j) = U(i,j);
            U(i,j) = 0.0;
        }
    }
}

} } }

#endif
