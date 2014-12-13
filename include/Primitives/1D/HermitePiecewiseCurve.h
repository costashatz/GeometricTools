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

#ifndef PIECEWISE_HERMITE_CURVE_H
#define PIECEWISE_HERMITE_CURVE_H

/**
* Includes
**/
#include <Math/Vector.h>
#include <Math/Matrix.h>
#include <Primitives/1D/HermiteCurve.h>
#include <Primitives/1D/PiecewiseCurve.h>
#include <vector>

using std::vector;


namespace LinearAlgebraTools {

using Math::Matrix;
using Math::Vector;

namespace Primitives {

/**
* HermiteCurve Class
* 1D (function) Piecewise Hermite Cubic curve
**/
class HermitePiecewiseCurve : public PiecewiseCurve
{
public:
    HermitePiecewiseCurve() {}

    void addDotPoint(const Vector<2> &point)
    {
        points.push_back(point[0]);
        dot_points.push_back(point);
        if(defined())
            calculateCoefficients();
    }

    void addDotPointInPlace(const Vector<2>& point, const unsigned int& i)
    {
        if(i>=points.size())
            return addDotPoint(point);
        int c = int(i)-1;
        if(c<0)
        {
            points.insert(points.begin(), point[0]);
            dot_points.insert(dot_points.begin(), point);
            curves.insert(curves.begin(), new HermiteCurve());
            ((HermiteCurve*)curves[0])->addDotPoint(dot_points[0]);
            ((HermiteCurve*)curves[0])->addDotPoint(dot_points[1]);
            return;
        }
        curves[c]->getPoints().clear();
        curves[c]->getDotPoints().clear();
        points.insert(points.begin()+i, point[0]);
        dot_points.insert(dot_points.begin()+i, point);
        ((HermiteCurve*)curves[c])->addDotPoint(dot_points[i-1]);
        ((HermiteCurve*)curves[c])->addDotPoint(dot_points[i]);
        curves.insert(curves.begin()+i, new HermiteCurve());
        ((HermiteCurve*)curves[i])->addDotPoint(dot_points[i]);
        ((HermiteCurve*)curves[i])->addDotPoint(dot_points[i+1]);
    }

protected:
    void calculateCoefficients()
    {
        if(!defined())
            return;
        for(int i=dot_points.size()-2;i<dot_points.size()-1;i++)
        {
            curves.push_back(new HermiteCurve());
            ((HermiteCurve*)curves[i])->addDotPoint(dot_points[i]);
            ((HermiteCurve*)curves[i])->addDotPoint(dot_points[i+1]);
        }
    }

    bool defined() const
    {
        return (points.size()>=2 && dot_points.size()>=2);
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
