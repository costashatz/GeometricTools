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

#ifndef GEOMETRIC_TOOLS_DISTANCES_POLYLINE_TO_POLYLINE_H
#define GEOMETRIC_TOOLS_DISTANCES_POLYLINE_TO_POLYLINE_H

/**
* Includes
**/
#include <geometric_tools/Primitives/LinearShapes.h>
#include <geometric_tools/Distances/LinearToPolyline.h>
#include <geometric_tools/Primitives/Polyline.h>

namespace GeometricTools {

using Primitives::Segment;
using Primitives::Polyline;

namespace Distances {

/**
* Computes Polyline to Polyline Distance Squared
* @param poly1
* @param poly2
**/
template<unsigned int N>
double distanceSq(const Polyline<N>& poly1, const Polyline<N>& poly2)
{
    double m = std::numeric_limits<double>::infinity();
    for(int i=0;i<poly1.vertices().size()-1;i++)
    {
        for(int j=0;j<poly2.vertices().size()-1;j++)
        {
            double tmp = distanceSq(Segment<N>(poly1.vertices()[i], poly1.vertices()[i+1]), Segment<N>(poly2.vertices()[j], poly2.vertices()[j+1]));
            if(tmp<m)
                m = tmp;
        }
    }
    return m;
}

/**
* Computes Polyline to Polyline Distance
* @param poly1
* @param poly2
**/
template<unsigned int N>
double distance(const Polyline<N>& poly1, const Polyline<N>& poly2)
{
    return sqrt(distanceSq(poly1, poly2));
}

} }

#endif
