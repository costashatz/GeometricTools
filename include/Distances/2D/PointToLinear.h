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

#ifndef POINT_TO_LINEAR_H
#define POINT_TO_LINEAR_H

/**
* Includes
**/
#include <Math/Vector.h>
#include <Primitives/LinearShapes.h>

namespace LinearAlgebraTools {

using Math::Vector;
using Primitives::Line;
using Primitives::Ray;
using Primitives::Segment;

namespace Distances {

/**
* Computes Point to Line Distance Squared
* @param point
* @param line
**/
double DistanceSq(const Vector<2>& point, const Line<2>& line)
{
    Vector<2> y_p = point-line.p();
    double k = (line.d()*y_p);
    return (y_p.lengthSq()-k*k);
}

double DistanceSq(const Line<2>& line, const Vector<2>& point)
{
    return DistanceSq(point,line);
}

/**
* Computes Point to Line Distance
* @param point
* @param line
**/
double Distance(const Vector<2>& point, const Line<2>& line)
{
    return sqrt(DistanceSq(point, line));
}

double Distance(const Line<2>& line, const Vector<2>& point)
{
    return Distance(point,line);
}

/**
* Computes Point to Ray Distance Squared
* @param point
* @param ray
**/
double DistanceSq(const Vector<2>& point, const Ray<2>& ray)
{
    Vector<2> toP = point-ray.p();
    if((ray.d()*toP)>0)
        return DistanceSq(point,Line<2>(ray.p(), ray.d()));
    return toP.lengthSq();
}

double DistanceSq(const Ray<2>& ray, const Vector<2>& point)
{
    return DistanceSq(point,ray);
}

/**
* Computes Point to Ray Distance
* @param point
* @param ray
**/
double Distance(const Vector<2>& point, const Ray<2>& ray)
{
    Vector<2> toP = point-ray.p();
    if((ray.d()*toP)>0)
        return Distance(point,Line<2>(ray.p(), ray.d()));
    return toP.length();
}

double Distance(const Ray<2>& ray, const Vector<2>& point)
{
    return Distance(point,ray);
}

/**
* Computes Point to Segment Distance Squared
* @param point
* @param seg
**/
double DistanceSq(const Vector<2>& point, const Segment<2>& seg)
{
    Vector<2> D = seg.d();
    Vector<2> toP = point-seg.p();
    double t = D*toP;
    if(t<=0)
    {
        return toP.lengthSq();
    }
    double DdD = D.lengthSq();
    if(t>=DdD)
    {
        Vector<2> toP1 = point-seg.P1();
        return toP1.lengthSq();
    }
    return DistanceSq(point, Line<2>(seg.P0(), seg.d()));
}

double DistanceSq(const Segment<2>& seg, const Vector<2>& point)
{
    return DistanceSq(point,seg);
}

/**
* Computes Point to Segment Distance
* @param point
* @param seg
**/
double Distance(const Vector<2>& point, const Segment<2>& seg)
{
    Vector<2> D = seg.d();
    Vector<2> toP = point-seg.p();
    double t = D*toP;
    if(t<=0)
    {
        return toP.length();
    }
    double DdD = D.lengthSq();
    if(t>=DdD)
    {
        Vector<2> toP1 = point-seg.P1();
        return toP1.length();
    }
    return Distance(point, Line<2>(seg.P0(), seg.d()));
}

double Distance(const Segment<2>& seg, const Vector<2>& point)
{
    return Distance(point,seg);
}

} }

#endif
