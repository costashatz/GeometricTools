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

#ifndef GEOMETRIC_TOOLS_PRIMITIVES_1D_POLYNOMIAL_CURVE_H
#define GEOMETRIC_TOOLS_PRIMITIVES_1D_POLYNOMIAL_CURVE_H

/**
* Includes
**/
#include <geometric_tools/Math/Vector.h>
#include <geometric_tools/Math/Matrix.h>
#include <geometric_tools/Math/LinearSystems/SolveLinear.h>
#include <geometric_tools/Primitives/1D/Curve.h>
#include <vector>
#include <map>

using std::vector;
using std::map;


namespace GeometricTools {

using Math::Matrix;
using Math::Vector;
using Math::LinearSystems::solveLinear;

namespace Primitives {

/**
* PolynomialCurve Class
* 1D (function) Polynomial curve defined by a series of control points (and slopes, curvatures) OR coefficients
**/
template<unsigned int N>
class PolynomialCurve : public Curve
{
protected:
    Vector<N+1> coefficients_;
    bool coeff_defined_;
    int n;
public:
    PolynomialCurve(): coeff_defined_(false), n(N) {}

    PolynomialCurve(const Vector<N+1>& coef): coefficients_(coef), coeff_defined_(true), n(N) {}

    virtual vector<double> coeff()
    {
        vector<double> tmp;
        tmp.insert(tmp.begin(), coefficients_.data(), coefficients_.data()+(n+1));
        return tmp;
    }

    virtual double getPoint(const double &u)
    {
        double s = 0.0;
        for(int i=0;i<=n;i++)
        {
            s += coefficients_[i]*std::pow(u, double(n-i));
        }
        return s;
    }

protected:
    void calculateCoefficients()
    {
        if(!defined())
            return;

        Matrix<N+1,N+1> A;
        Vector<N+1> y;
        map<double,double> us;
        for(int i=0;i<points_.size();i++)
        {
            double d = double(points_.size()-1);
            if(points_.size()<=1)
                d = 1.0;
            us[points_[i]] = double(i)/d;
        }
        int ddot = ddot_points_.size();
        int dot = dot_points_.size();
        for(int i=0;i<ddot;i++)
        {
            y[i] = ddot_points_[i][1];
            for(int j=0;j<=n-2;j++)
            {
                A(i,j) = (n-j)*(n-j-1)*std::pow(us[ddot_points_[i][0]], n-j-2);
            }
        }
        for(int i=ddot;i<ddot+dot;i++)
        {
            y[i] = dot_points_[i-ddot][1];
            for(int j=0;j<=n-1;j++)
            {
                A(i,j) = (n-j)*std::pow(us[dot_points_[i-ddot][0]], n-j-1);
            }
        }
        for(unsigned int i=dot+ddot;i<=n;i++)
        {
            y[i] = points_[i-(dot+ddot)];
            for(unsigned int j=0;j<=n;j++)
            {
                A(i,j) = std::pow(us[points_[i-(dot+ddot)]], n-j);
            }
        }
        coefficients_ = solveLinear(A, y);
    }

    bool defined() const
    {
        if(coeff_defined_)
            return true;
        if (points_.size() + dot_points_.size() + ddot_points_.size() <= N)
            return false;
        for(int i=0;i<=n;i++)
        {
            int d = dot_points_.size();
            int dd = ddot_points_.size();
            if(i>=n)
                d = 0;
            if(i>=(n-1))
                dd = 0;
            if((points_.size()+d+dd)==0)
                return false;
        }
        return true;
    }

    bool canAddPoint(const double& point)
    {
        return true;
    }

    bool canAddDotPoint(const Vector<2>& point)
    {
        return int(points_.size()+dot_points_.size()+ddot_points_.size()) <= (n-1);
    }

    bool canAddDDotPoint(const Vector<2>& point)
    {
        return int(points_.size()+dot_points_.size()+ddot_points_.size()) <= (n-1);
    }
};

/**
* Typedefs for frequently used types
**/
typedef PolynomialCurve<0> ConstantPolynomialCurve;

typedef PolynomialCurve<1> LinearPolynomialCurve;

typedef PolynomialCurve<2> QuadraticPolynomialCurve;

typedef PolynomialCurve<3> CubicPolynomialCurve;

} }

#endif
