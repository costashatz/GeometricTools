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

#ifndef GOLDEN_SEARCH_MINIMIZATION_H
#define GOLDEN_SEARCH_MINIMIZATION_H

/**
* Includes
**/
#include <tuple>
#include <cmath>

namespace GeometricTools { namespace Math { namespace NumericalOptimization {

inline std::tuple<double,double,double> getBracket(double(*f)(double), const double& x1, const double& x2)
{
    double X1 = (x1<x2)? x1 : x2;
    double X2 = (x2>x1)? x2 : x1;
    double d = sqrt(X2-X1);
    bool cond = f(X1)<f(X2);
    double X3 = (cond)? X1-d : X2+d;
    if(f(X3)>std::min(f(X1), f(X2)))
    {
        if(cond)
            return std::make_tuple(X3, X1, X2);
        return std::make_tuple(X1, X2, X3);
    }
    X2 = (f(X1)<f(X2))? X1 : X2;
    X1 = X3;
    return getBracket(f, X1, X2);
}

inline std::tuple<double,double> goldenSearchMinimize(double(*f)(double), const double& a, const double& b, const double& c, const double& tau_ = 1e-10)
{
    double tau = sqrt(tau_);
    static double w_ = 2-(1+sqrt(5))/2;
    double x;
    if((c-b)>(b-a))
        x = b+w_*(c-b);
    else
        x = b-w_*(b-a);
    if (std::abs(c-a) < tau*(std::abs(b)+std::abs(x)))
        return std::make_tuple((c+a)/2.0, f((c+a)/2.0));
    if(f(x)<f(b))
    {
        if((c-b) > (b-a))
            return goldenSearchMinimize(f, b, x, c, tau);
        else
            return goldenSearchMinimize(f, a, x, b, tau);
    }
    else
    {
        if((c-b) > (b-a))
            return goldenSearchMinimize(f, a, b, x, tau);
        else
            return goldenSearchMinimize(f, x, b, c, tau);
    }
}

} } }

#endif
