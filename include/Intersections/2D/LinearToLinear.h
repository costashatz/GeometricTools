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

#ifndef GEOMETRIC_TOOLS_INTERSECTIONS_2D_LINEAR_TO_LINEAR_H
#define GEOMETRIC_TOOLS_INTERSECTIONS_2D_LINEAR_TO_LINEAR_H

/**
* Includes
**/
#include <Math/Vector.h>
#include <Primitives/LinearShapes.h>
#include <Intersections/IntersectionInfo.h>
#include <limits>

namespace GeometricTools {

using Math::Vector;
using Primitives::Line;
using Primitives::Ray;
using Primitives::Segment;

namespace Intersections {

inline Intersection2DInfo* intersect(const Segment<2>& seg1, const Segment<2>& seg2)
{
    Intersection2DInfo* info = new Intersection2DInfo;
    double epsilon = std::numeric_limits<double>::epsilon();
    Vector<2> u = seg1.d();
    Vector<2> v = seg2.d();
    Vector<2> w = seg1.P0()-seg2.P0();
    double D = u[0]*v[1]-u[1]*v[0];
    double tmp1 = u[0]*w[1]-u[1]*w[0], tmp2 = v[0]*w[1]-v[1]*w[0];
    if(std::abs(D)<epsilon)
    {
        if(std::abs(tmp1)>epsilon || std::abs(tmp2)>epsilon)
            return nullptr;
        double du = u*u;
        double dv = v*v;
        if(du<epsilon && dv<epsilon)
        {
            if(std::abs(du-dv)>epsilon)
                return nullptr;
            info->point = seg1.P0();
            return info;
        }
        if(du<epsilon)
        {
            //TODO: test for inclusion of seg1.P0 in the seg2
            info->point = seg1.P0();
            return info;
        }
        if(dv<epsilon)
        {
            //TODO: test for inclusion of seg2.P0 in the seg1
            info->point = seg2.P0();
            return info;
        }
        double t0, t1;
        Vector<2> w2 = seg1.P1()-seg2.P0();
        if(std::abs(v[0])>epsilon)
        {
            t0 = w[0]/v[0];
            t1 = w2[0]/v[0];
        }
        else
        {
            t0 = w[1]/v[1];
            t1 = w2[1]/v[1];
        }
        if(t0>t1)
        {
            double t = t0;
            t0 = t1;
            t1 = t0;
        }
        if(t0>1 || t1 < 0)
            return nullptr;
        t0 = (t0<0)? 0 : t0;
        t1 = (t1>1)? 1 : t1;
        if((t1-t0)<epsilon)
        {
            info->point = seg2.P0()+t0*v;
            return info;
        }
        info->point = seg2.P0()+t0*v;
        info->delta = (seg2.P0()+t1*v)-info->point;
        return info;
    }
    double sI = tmp1/D;
    if(sI<0.0 || sI>1.0)
        return nullptr;
    double tI = tmp2/D;
    if(tI<0.0 || tI>1.0)
        return nullptr;

    info->point = seg1.P0()+sI*u;
    return info;
}

} }


#endif
