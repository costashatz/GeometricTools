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

#ifndef GEOMETRIC_TOOLS_PRIMITIVES_1D_HERMITE_CURVE_H
#define GEOMETRIC_TOOLS_PRIMITIVES_1D_HERMITE_CURVE_H

/**
* Includes
**/
#include <Math/Vector.h>
#include <Math/Matrix.h>
#include <Math/LinearSystems/SolveLU.h>
#include <Primitives/1D/Curve.h>
#include <vector>

using std::vector;


namespace GeometricTools {

using Math::Matrix;
using Math::Vector;
using Math::LinearSystems::solveLU;

namespace Primitives {

/**
* HermiteCurve Class
* 1D (function) Hermite Cubic curve
**/
class HermiteCurve : public Curve
{
protected:
    Vector<4> coefficients_;
    bool coeff_defined_;
public:
    HermiteCurve(): coeff_defined_(false) {}

    HermiteCurve(const Vector<4>& coef): coefficients_(coef), coeff_defined_(true) {}

    virtual vector<double> coeff()
    {
        vector<double> tmp;
        tmp.insert(tmp.begin(), coefficients_.data(), coefficients_.data()+4);
        return tmp;
    }

    virtual double getPoint(const double &u)
    {
        double s = 0.0;
        for(int i=0;i<=3;i++)
        {
            s += coefficients_[i]*std::pow(u, double(3-i));
        }
        return s;
    }

protected:
    void calculateCoefficients()
    {
        if(!defined())
            return;
        double p0 = points_[0], pd0 = dot_points_[0][1];
        double p1 = points_[1], pd1 = dot_points_[1][1];
        Matrix<4,4> A{0.0, 0.0, 0.0, 1.0,
                      1.0, 1.0, 1.0, 1.0,
                      0.0, 0.0, 1.0, 0.0,
                      3.0, 2.0, 1.0, 0.0};
        if(p0!=p1)
        {
            pd0 /= (p1-p0);
            pd1 /= (p1-p0);
        }
        Vector<4> y{p0, p1, pd0, pd1};
        coefficients_ = solveLU(A, y);
    }

    bool defined() const
    {
        if(coeff_defined_)
            return true;
        return (points_.size()==2 && dot_points_.size()==2);
    }

    bool canAddPoint(const double& point)
    {
        return false;
    }

    bool canAddDotPoint(const Vector<2>& point)
    {
        return true;
    }

    bool canAddDDotPoint(const Vector<2>& point)
    {
        return false;
    }
};

} }

#endif
