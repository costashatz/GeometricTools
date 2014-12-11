#ifndef QR_DECOMPOSITION
#define QR_DECOMPOSITION

/**
* Includes
**/
#include <Math/Matrix.h>

extern "C" int dgeqrf_(int *m, int *n, double *a, int *
                       lda, double *tau, double *work, int *lwork, int *info);

extern "C" int dorgqr_(int *m, int *n, int *k, double *
                       a, int *lda, double *tau, double *work, int *lwork,
                       int *info);

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
    int dim = int(D);
    double* work = new double[D];
    dgeqrf_(&dim, &dim, R.data(), &dim, tau, work, &dim, info);
    Q = R;
    for(unsigned int i=0;i<D;i++)
    {
        for(unsigned int j=0;j<i;j++)
        {
            R(i,j) = 0.0;
        }
    }
    dorgqr_(&dim, &dim, &dim, Q.data(), &dim, tau, work, &dim, &info);
}

} } }

#endif
