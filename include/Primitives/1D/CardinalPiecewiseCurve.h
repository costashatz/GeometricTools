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

#ifndef PIECEWISE_CARDINAL_CURVE_H
#define PIECEWISE_CARDINAL_CURVE_H

/**
* Includes
**/
#include <Math/Vector.h>
#include <Math/Matrix.h>
#include <Primitives/1D/CardinalCurve.h>
#include <Primitives/1D/PiecewiseCurve.h>
#include <vector>

using std::vector;


namespace GeometricTools {

using Math::Matrix;
using Math::Vector;

namespace Primitives {

/**
* CardinalPiecewiseCurve Class
* 1D (function) Piecewise Cardinal Cubic curve
**/
class CardinalPiecewiseCurve : public PiecewiseCurve
{
protected:
    double t_;
public:
    CardinalPiecewiseCurve(const double& t = 0.0): t_(t) {}

    double& t()
    {
        return t_;
    }

    void setT(const double& t)
    {
        t_ = t;
    }

    void addPoint(const double& point)
    {
        points.push_back(point);
        if(defined())
            calculateCoefficients();
    }

protected:
    void calculateCoefficients()
    {
        if(!defined())
            return;
        curves.push_back(new CardinalCurve(t_));
        for(int i=points.size()-4;i<points.size();i++)
        {
            ((CardinalCurve*)curves[curves.size()-1])->addPoint(points[i]);
        }
    }

    bool defined() const
    {
        return (points.size()>=4);
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
