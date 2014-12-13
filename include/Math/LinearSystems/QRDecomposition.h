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

namespace GeometricTools {

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
