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

#ifndef GEOMETRIC_TOOLS_PRIMITIVES_2D_CIRCLE_H
#define GEOMETRIC_TOOLS_PRIMITIVES_2D_CIRCLE_H

/**
* Includes
**/
#include <geometric_tools/Misc/Helper.h>
#include <geometric_tools/Math/Matrix.h>
#include <geometric_tools/Primitives/2D/QuadraticCurve.h>


namespace GeometricTools {

using Math::Vector;

namespace Primitives {

/**
* Circle Class
* Circle defined by center and radius
**/
class Circle
{
protected:
    Vector<2> center_;
    double radius_;
public:
    Circle() {}

    Circle(const Vector<2>& center, const double& radius): center_(center), radius_(radius) {}

    Vector<2> center() { return center_; }

    double radius() { return radius_; }

    double area() { return Helper::Pi*radius_*radius_; }

    QuadraticCurve asQuadratic() { return QuadraticCurve({}, -2.0*center_, center_*center_-radius_*radius_); }
};

} }

#endif
