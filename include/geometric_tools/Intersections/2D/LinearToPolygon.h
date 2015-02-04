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

#ifndef GEOMETRIC_TOOLS_INTERSECTIONS_2D_LINEAR_TO_POLYGON_H
#define GEOMETRIC_TOOLS_INTERSECTIONS_2D_LINEAR_TO_POLYGON_H

/**
* Includes
**/
#include <geometric_tools/Math/Vector.h>
#include <geometric_tools/Primitives/LinearShapes.h>
#include <geometric_tools/Primitives/2D/Polygon.h>
#include <geometric_tools/Primitives/Polyline.h>
#include <geometric_tools/Intersections/IntersectionInfo.h>
#include <limits>

namespace GeometricTools {

using Math::Vector;
using Primitives::Line;
using Primitives::Ray;
using Primitives::Segment;
using Primitives::Polygon;
using Primitives::Polyline;

namespace Intersections {

inline Intersection2DInfo* intersect(const Segment<2>& seg, const Polygon& poly)
{
    Intersection2DInfo* info = new Intersection2DInfo;
    if(seg.P0()==seg.P1())
    {
        //TODO: test for inclusion of seg.P0 in the poly
        info->point = seg.P0();
        info->delta = seg.d();
        return info;
    }
    double tE = 0.0, tL = 1.0;
    double t, N, D;
    int n = poly.vertices().size();
    Vector<2> dS = seg.d();
    for(int i=0;i<n;i++)
    {
        int i_p = (i+1)%n;
        Vector<2> e = poly.vertices()[i_p]-poly.vertices()[i];
        Vector<2> tmp = seg.P0()-poly.vertices()[i];
        N = e[0]*tmp[1]-e[1]*tmp[0];
        D = -(e[0]*dS[1]-e[1]*dS[0]);
        if(std::abs(D)<std::numeric_limits<double>::epsilon())
        {
            if(N<0)
                return nullptr;
            else
                continue;
        }
        t = N/D;
        if(D<0)
        {
            if(t>tE)
            {
                tE = t;
                if(tE>tL)
                    return nullptr;
            }
        }
        else
        {
            if(t<tL)
            {
                tL = t;
                if(tL<tE)
                    return nullptr;
            }
        }
    }
    info->point = seg.P0()+tE*dS;
    info->delta = (seg.P0()+tL*dS)-info->point;
    return info;
}

inline Intersection2DInfo* intersect(const Polygon& poly, const Segment<2>& seg)
{
    return intersect(seg, poly);
}

inline Intersection2DInfo* intersect(const Polyline<2>& line, const Polygon& poly)
{
    Intersection2DInfo* firstInfo = nullptr;
    Intersection2DInfo* lastInfo = nullptr;
    for(int i=0;i<line.vertices().size()-1;i++)
    {
        Intersection2DInfo* info = intersect(Segment<2>(line.vertices()[i], line.vertices()[i+1]), poly);
        if(info!=nullptr && firstInfo!=nullptr)
            firstInfo = info;
        else if(info!=nullptr)
            lastInfo = info;
    }
    if(firstInfo==nullptr && lastInfo==nullptr)
        return nullptr;
    if(lastInfo==nullptr)
        return firstInfo;
    Vector<2> point = lastInfo->point+lastInfo->delta;
    firstInfo->delta = (point-firstInfo->point);
    return firstInfo;
}

} }

#endif
