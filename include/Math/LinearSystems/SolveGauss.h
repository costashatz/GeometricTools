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

#ifndef SOLVE_GAUSS_H
#define SOLVE_GAUSS_H

/**
* Includes
**/
#include <limits>

namespace LinearAlgebraTools {

namespace Math {

template<unsigned int ROWS, unsigned int COLS>
class Matrix;

template<unsigned int N>
class Vector;

namespace LinearSystems {

/**
* Solve Linear System using Gauss Elimination
* @param A - parameter Matrix
* @param B - constant Vector
* @return Vector<D> - solution
**/
template<unsigned int D>
Vector<D> solveGauss(const Matrix<D,D>& A, const Vector<D>& B)
{
    Matrix<D,D> a = A;
    Vector<D> b = B;
    Vector<D> res;
    for(unsigned int i=0;i<D-1;i++)
    {
        int m = 0;
        double temp1 = 0.0;
        for(unsigned int j=i;j<D;j++)
        {
            double temp2 = 0.0;
            for(unsigned int k=i+1;k<D;k++)
                if(abs(a(j,k))>temp2)
                    temp2 = abs(a(j,k));
            if(abs(a(j,i))>(temp1*temp2))
            {
                temp1 = abs(a(j,i))/temp2;
                m = j;
            }
        }
        if(m!=i)
        {
            a.swapRows(i,m);
            double temp = b(i);
            b(i) = b(m);
            b(m) = temp;
        }

        if(abs(a(i,i))<std::numeric_limits<double>::epsilon())
            return res;

        for(unsigned int k=i+1;k<D;k++)
        {
            double r = a(k,i)/a(i,i);
            for(unsigned int c=0;c<D;c++)
            {
                a(k,c) -= r*a(i,c);
            }
            b(k) -= r*b(i);
        }
    }
    for(int k=D-1;k>=0;k--)
    {
        double S = 0.0;
        for(unsigned int j=k+1;j<D;j++)
        {
            S += a(k,j)*res(j);
        }
        res(k) = 1.0/a(k,k)*(b(k)-S);
    }
    return res;
}

} } }

#endif
