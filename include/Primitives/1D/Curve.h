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

#ifndef GEOMETRIC_TOOLS_PRIMITIVES_1D_CURVE_H
#define GEOMETRIC_TOOLS_PRIMITIVES_1D_CURVE_H

/**
* Includes
**/
#include <Math/Vector.h>
#include <Misc/Helper.h>
#include <vector>

using std::vector;

namespace GeometricTools {

using Math::Vector;
using namespace Helper;

namespace Primitives {

/**
* Curve Class
* Generic 1D (function) curve - base class for all parametric curves
**/
class Curve
{
protected:
    vector<double> points;
    vector<Vector<2> > dot_points;
    vector<Vector<2> > ddot_points;
protected:
    virtual void calculateCoefficients() = 0;

    virtual bool defined() const = 0;

    virtual bool canAddPoint(const double& point) = 0;

    virtual bool canAddDotPoint(const Vector<2>& point) = 0;

    virtual bool canAddDDotPoint(const Vector<2>& point) = 0;
public:
    void addPoint(const double& point)
    {
        if(defined() || !canAddPoint(point))
            return;
        points.push_back(point);
        if(defined())
            calculateCoefficients();
    }

    void addPointInPlace(const double& point, const unsigned int& i)
    {
        //Override to give desired functionality
        addPoint(point);
    }

    void addDotPoint(const Vector<2>& point)
    {
        if(defined() || !canAddDotPoint(point))
            return;
        points.push_back(point[0]);
        dot_points.push_back(point);
        if(defined())
            calculateCoefficients();
    }

    void addDotPointInPlace(const Vector<2>& point, const unsigned int& i)
    {
        //Override to give desired functionality
        addDotPoint(point);
    }

    void addDDotPoint(const Vector<2>& point)
    {
        if(defined() || !canAddDDotPoint(point))
            return;
        points.push_back(point[0]);
        ddot_points.push_back(point);
        if(defined())
            calculateCoefficients();
    }

    void addDDotPointInPlace(const Vector<2>& point, const unsigned int& i)
    {
        //Override to give desired functionality
        addDDotPoint(point);
    }

    virtual double getPoint(const double& u) = 0;

    virtual vector<double> coeff() = 0;

    vector<double>& getPoints() { return points; }
    vector<Vector<2> >& getDotPoints() { return dot_points; }
    vector<Vector<2> >& getDDotPoints() { return ddot_points; }
};

} }

#endif
