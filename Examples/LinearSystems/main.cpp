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

#include <iostream>
#include <Math/LinearSystems/SolveGauss.h>
#include <Math/LinearSystems/SolveLU.h>
#include <Math/LinearSystems/SolveLinear.h>
using namespace std;

using namespace GeometricTools::Math;
using namespace GeometricTools::Math::LinearSystems;

//SIMPLE EXAMPLE

int main(int argc, char *argv[])
{
    Matrix<3,3> A = Matrix<3,3>(3.,-6.,-3,2.,0.,6.,-4.,7.,4.);
//    Matrix<4,4> A = Matrix<4,4>(5., 7., 6., 5.,
//                                7., 10., 8., 7,
//                                6., 8., 10., 9.,
//                                5., 7., 9., 10.);
    Vector<3> B = Vector<3>(-3.,-22.,3.);
    Vector<3> x = solveLU(A,B);
    cout<<"Solving System: \n";
    cout<<A<<endl;
    cout<<" = \n";
    cout<<B<<endl;
    cout<<"\nSOLUTION: ";
    cout<<x<<endl;
    cout<<"------------------------------------\n\n";
    Matrix<4,4> a = Matrix<4,4>(2,1,1,0,4,3,3,1,8,7,9,5,6,7,9,8);
    Vector<4> b = Vector<4>(1.,2.,3,4);
    cout<<"SOLUTION with LU: ";
    Matrix<4,4> a1 = a;
    Vector<4> y = solveLU(a,b);
    cout<<y<<endl;
    Vector<4> y2 = solveGauss(a,b);
    cout<<"SOLUTION with Gauss: ";
    cout<<y2<<endl;
    cout<<"------------------------------------\n\n";
    Vector<4> y3 = solveLinear(a,b);
    cout<<"SOLUTION with Linear: ";
    cout<<y3<<endl;
    cout<<"------------------------------------\n\n";
    return 0;
}
