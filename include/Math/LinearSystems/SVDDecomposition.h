#ifndef SVD_DECOMPOSITION
#define SVD_DECOMPOSITION

/**
* Includes
**/
#include <Math/Matrix.h>

extern "C" int dgesvd_(char *jobu, char *jobvt, int *m, int *n,
                       double *a, int *lda, double *s, double *u, int *
                       ldu, double *vt, int *ldvt, double *work, int *lwork,
                       int *info);

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
    double* work = new double[D];
    int dim = int(D), info, lwork = std::max(4*D, 5*D);
    char job = 'S';
    //TODO: DOES NOT WORK
    dgesvd_(&job, &job, &dim, &dim, tmp.data(), &dim, S.data(), U.data(), &dim, V.data(), &dim, work, &lwork, &info);
    U = U.transpose();
    S = S.transpose();
    for(unsigned int i=1;i<D;i++) {
        S(i,i) = S(0,i);
        S(0,i) = 0.0;
    }
}

} } }

#endif
