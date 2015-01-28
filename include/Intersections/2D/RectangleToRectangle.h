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

#ifndef GEOMETRIC_TOOLS_INTERSECTIONS_2D_RECTANGLE_TO_RECTANGLE_H
#define GEOMETRIC_TOOLS_INTERSECTIONS_2D_RECTANGLE_TO_RECTANGLE_H

/**
* Includes
**/
#include <Math/Vector.h>
#include <Primitives/2D/Rectangle.h>
#include <Intersections/IntersectionInfo.h>
#include <limits>

namespace GeometricTools {

using Math::Vector;
using Primitives::Rectangle;

namespace Intersections {

inline Intersection2DInfo* intersect(const Rectangle& r1, const Rectangle& r2)
{
    Intersection2DInfo* info = new Intersection2DInfo;

    Vector<2> center1 = r1.center();

    Vector<2> center2 = r2.center();

    Vector<2> half1 = r1.half();
    Vector<2> half2 = r2.half();

    double half1x = half1[0];
    double half1y = half1[1];

    double half2x = half2[0];
    double half2y = half2[1];

    Vector<2> dC = center2-center1;
    double px = half1x+half2x-std::abs(dC[0]);
    if(px < 0)
        return nullptr;
    double py = half1y+half2y-std::abs(dC[1]);
    if(py < 0)
        return nullptr;
    Vector<2> p;
    Vector<2> n;
    if(px<py)
    {
        int sx = std::abs(dC[0])/dC[0];
        p = Vector<2>(center1[0]+half1x*sx, center2[1]);
        n = Vector<2>(sx*px, 0.0);
    }
    else
    {
        int sy = std::abs(dC[1])/dC[1];
        p = Vector<2>(center2[0], center1[1]+half1y*sy);
        n = Vector<2>(0.0, sy*py);
    }
    info->point = p;
    info->delta = n;
    return info;
}

} }

#endif
