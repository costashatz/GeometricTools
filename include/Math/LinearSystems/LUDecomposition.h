#ifndef LU_DECOMPOSITION
#define LU_DECOMPOSITION

/**
* Includes
**/
#include <Math/Matrix.h>
#include <lapacke.h>

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
    LAPACKE_dgetrf(LAPACK_ROW_MAJOR, D, D, U.values, D, ipiv);
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
