#ifndef QR_DECOMPOSITION
#define QR_DECOMPOSITION

/**
* Includes
**/
#include <Math/Matrix.h>
#include <lapacke.h>

namespace LinearAlgebraTools {

namespace Math {

namespace LinearSystems {

/**
* Perform QR Decomposition
* @param a - Matrix to decompose
* @param Q - Q Matrix --
* @param R - R Matrix  |
*					   --> Pass by reference (returns)
**/
template<unsigned int D>
void QRDecomposition(const Matrix<D,D>& a, Matrix<D,D>& Q, Matrix<D,D>& R)
{
    Q.identity();
    R = a;
    double* tau = new double[D];
    LAPACKE_dgeqrf(LAPACK_ROW_MAJOR, D, D, R.data(), D, tau);
    Q = R;
    for(unsigned int i=0;i<D;i++)
    {
        for(unsigned int j=0;j<i;j++)
        {
            R(i,j) = 0.0;
        }
    }
    LAPACKE_dorgqr(LAPACK_ROW_MAJOR, D, D, D, Q.data(), D, tau);
}

} } }

#endif
