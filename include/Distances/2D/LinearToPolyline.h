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

#ifndef LINEAR_TO_POLYLINE_H
#define LINEAR_TO_POLYLINE_H

/**
* Includes
**/
#include <Math/Vector.h>
#include <Primitives/LinearShapes.h>
#include <Distances/2D/PointToPolyline.h>
#include <Distances/2D/LinearToLinear.h>
#include <Primitives/2D/Polyline.h>

namespace LinearAlgebraTools {

using Math::Vector;
using Primitives::Line;
using Primitives::Ray;
using Primitives::Segment;
using Primitives::Polyline;

namespace Distances {


/**
* Computes Segment to Polyline Distance Squared
* @param segment
* @param polyline
**/
double DistanceSq(const Segment<2>& seg, const Polyline& polyline)
{
    double m = DistanceSq(seg, Segment<2>(polyline.vertices()[0], polyline.vertices()[1]));
    for(unsigned int i=1;i<polyline.vertices().size()-1;i++)
    {
        double tmp = DistanceSq(seg, Segment<2>(polyline.vertices()[i], polyline.vertices()[i+1]));
        if(tmp<m)
            m = tmp;
    }
    return m;
}

double DistanceSq(const Polyline& polyline, const Segment<2>& seg)
{
    return DistanceSq(seg, polyline);
}


/**
* Computes Segment to Polyline Distance
* @param segment
* @param polyline
**/
double Distance(const Segment<2>& seg, const Polyline& polyline)
{
    double m = Distance(seg, Segment<2>(polyline.vertices()[0], polyline.vertices()[1]));
    for(unsigned int i=1;i<polyline.vertices().size()-1;i++)
    {
        double tmp = Distance(seg, Segment<2>(polyline.vertices()[i], polyline.vertices()[i+1]));
        if(tmp<m)
            m = tmp;
    }
    return m;
}

double Distance(const Polyline& polyline, const Segment<2>& seg)
{
    return Distance(seg, polyline);
}

} }

#endif
