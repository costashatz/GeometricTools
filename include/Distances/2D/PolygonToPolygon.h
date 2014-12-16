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

#ifndef GEOMETRIC_TOOLS_DISTANCES_2D_POLYGON_TO_POLYGON_H
#define GEOMETRIC_TOOLS_DISTANCES_2D_POLYGON_TO_POLYGON_H

/**
* Includes
**/
#include <Distances/PolylineToPolyline.h>
#include <Primitives/Polyline.h>
#include <Primitives/2D/Polygon.h>

namespace GeometricTools {

using Primitives::Polygon;
using Primitives::Polyline;

namespace Distances {

/**
* Computes Polygon to Polygon Distance Squared
* @param poly1
* @param poly2
**/
inline double distanceSq(const Polygon& poly1, const Polygon& poly2)
{
    double m = std::numeric_limits<double>::infinity();
    for(int i=0;i<poly1.vertices().size();i++)
    {
        int i_p = (i+1)%poly1.vertices().size();
        for(int j=0;j<poly2.vertices().size();j++)
        {
            int j_p = (j+1)%poly2.vertices().size();
            double tmp = distanceSq(Segment<2>(poly1.vertices()[i], poly1.vertices()[i_p]), Segment<2>(poly2.vertices()[j], poly2.vertices()[j_p]));
            if(tmp<m)
                m = tmp;
        }
    }
    return m;
}

/**
* Computes Polygon to Polygon Distance
* @param poly1
* @param poly2
**/
inline double distance(const Polygon& poly1, const Polygon& poly2)
{
    return sqrt(distanceSq(poly1, poly2));
}

} }

#endif
