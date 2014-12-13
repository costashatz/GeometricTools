/**
Copyright (c) 2014, Konstantinos Chatzilygeroudis
All rights reserved.
Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer
    in the documentation and/or other materials provided with the distribution.
3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived
    from this software without specific prior written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING,
BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
**/

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

namespace GeometricTools {

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
