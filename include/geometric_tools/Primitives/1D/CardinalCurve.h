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

#ifndef GEOMETRIC_TOOLS_PRIMITIVES_1D_CARDINAL_CURVE_H
#define GEOMETRIC_TOOLS_PRIMITIVES_1D_CARDINAL_CURVE_H

/**
* Includes
**/
#include <geometric_tools/Math/Vector.h>
#include <geometric_tools/Math/Matrix.h>
#include <geometric_tools/Math/LinearSystems/SolveLU.h>
#include <geometric_tools/Primitives/1D/Curve.h>
#include <vector>

using std::vector;


namespace GeometricTools {

using Math::Matrix;
using Math::Vector;
using Math::LinearSystems::solveLU;

namespace Primitives {

/**
* CardinalCurve Class
* 1D (function) Cardinal Cubic curve
**/
class CardinalCurve : public Curve
{
protected:
    Vector<4> coefficients_;
    bool coeff_defined_;
    double t_;
public:
    CardinalCurve(const double& t = 0.0): t_(t), coeff_defined_(false) {}

    CardinalCurve(const Vector<4>& coef, const double& t = 0.0): t_(t), coefficients_(coef), coeff_defined_(true) {}

    double& t()
    {
        return t_;
    }

    void setT(const double& t)
    {
        t_ = t;
    }

    virtual const double getPoint(const double &u) const
    {
        double s = 0.0;
        for(int i=0;i<=3;i++)
        {
            s += coefficients_[i]*std::pow(u, double(3-i));
        }
        return s;
    }

    virtual vector<double> coeff()
    {
        vector<double> tmp;
        tmp.insert(tmp.begin(), coefficients_.data(), coefficients_.data()+4);
        return tmp;
    }

protected:
    void calculateCoefficients()
    {
        if(!defined())
            return;
        double p0 = points_[1], pd0 = 0.5*(1.0-t_)*(points_[2]-points_[0]);
        double p1 = points_[2], pd1 = 0.5*(1.0-t_)*(points_[3]-points_[1]);
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
        return (points_.size()==4);
    }

    bool canAddPoint(const double& point)
    {
        return true;
    }

    bool canAddDotPoint(const Vector<2>& point)
    {
        return false;
    }

    bool canAddDDotPoint(const Vector<2>& point)
    {
        return false;
    }
};

} }

#endif
