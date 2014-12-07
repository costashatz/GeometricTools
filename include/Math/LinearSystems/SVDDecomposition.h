#ifndef SVD_DECOMPOSITION
#define SVD_DECOMPOSITION

/**
* Includes
**/
#include <Math/Matrix.h>
#include <lapacke.h>

namespace LinearAlgebraTools {

namespace Math {

namespace LinearSystems {

/**
* Perform SVD Decomposition
* @param a - Matrix to decompose
* @param U - U Matrix --
* @param S - S Matrix
* @param V - V Matrix  |
*					   --> Pass by reference (returns)
**/
template<unsigned int D>
void SVDDecomposition(const Matrix<D,D>& a, Matrix<D,D>& U, Matrix<D,D>& S, Matrix<D,D>& V)
{
    Matrix<D,D> tmp = a;
    U.identity();
    S.identity();
    V.identity();
    double* sup = new double[D];
    LAPACKE_dgesvd(LAPACK_ROW_MAJOR, 'S', 'S', D, D, tmp.values, D, S.values, U.values, D, V.values, D, sup);
    V = V.transpose();
    for(unsigned int i=1;i<D;i++) {
        S(i,i) = S(0,i);
        S(0,i) = 0.0;
    }
}

} } }

#endif
