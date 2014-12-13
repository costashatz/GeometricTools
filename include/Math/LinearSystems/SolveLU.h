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

#ifndef SOLVE_LU_H
#define SOLVE_LU_H

/**
* Includes
**/
#include <Math/Matrix.h>

extern "C" int dgetrf_(int *m, int *n, double *a, int * lda, int *ipiv, int *info);
extern "C" int dgetrs_(char *trans, int *n, int *nrhs,
                       double *a, int *lda, int *ipiv, double *b, int * ldb, int *info);

namespace GeometricTools {

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
