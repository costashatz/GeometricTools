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
