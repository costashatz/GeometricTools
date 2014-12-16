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
#include <Math/Vector.h>
#include <Primitives/LinearShapes.h>
#include <Primitives/2D/Polygon.h>
#include <limits>
#include <tuple>
#include <vector>

namespace GeometricTools {

using Math::Vector;
using Primitives::Line;
using Primitives::Ray;
using Primitives::Segment;
using Primitives::Polygon;

namespace Intersections {

using std::vector;

inline std::tuple<int, vector<Vector<2> > > intersect(const Segment<2>& seg, const Polygon& poly)
{
    vector<Vector<2> > points;
    if(seg.P0()==seg.P1())
    {
        //TODO: test for inclusion of seg.P0 in the poly
        points.push_back(seg.P0());
        points.push_back(seg.P1());
        return std::make_tuple(2, points);
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
                return std::make_tuple(0, points);
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
                    return std::make_tuple(0, points);
            }
        }
        else
        {
            if(t<tL)
            {
                tL = t;
                if(tL<tE)
                    return std::make_tuple(0, points);
            }
        }
    }
    points.push_back(seg.P0()+tE*dS);
    points.push_back(seg.P0()+tL*dS);
    return std::make_tuple(2, points);
}

inline std::tuple<int, vector<Vector<2> > > intersect(const Polygon& poly, const Segment<2>& seg)
{
    return intersect(seg, poly);
}

} }

#endif
