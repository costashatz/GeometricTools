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

#ifndef GEOMETRIC_TOOLS_DISTANCES_POINT_TO_LINEAR_H
#define GEOMETRIC_TOOLS_DISTANCES_POINT_TO_LINEAR_H

/**
* Includes
**/
#include <geometric_tools/Math/Vector.h>
#include <geometric_tools/Primitives/LinearShapes.h>

namespace GeometricTools {

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
template<unsigned int N>
double distanceSq(const Vector<N>& point, const Line<N>& line)
{
    Vector<N> w = point-line.p();
    Vector<N> u = line.d().normalized();
    return (w-(w*u)*u).lengthSq();
}

template<unsigned int N>
double distanceSq(const Line<N>& line, const Vector<N>& point)
{
    return distanceSq(point,line);
}

/**
* Computes Point to Line Distance
* @param point
* @param line
**/
template<unsigned int N>
double distance(const Vector<N>& point, const Line<N>& line)
{
    Vector<N> w = point-line.p();
    Vector<N> u = line.d().normalized();
    return (w-(w*u)*u).length();
}

template<unsigned int N>
double distance(const Line<N>& line, const Vector<N>& point)
{
    return distance(point,line);
}

/**
* Computes Point to Ray Distance Squared
* @param point
* @param ray
**/
template<unsigned int N>
double distanceSq(const Vector<N>& point, const Ray<N>& ray)
{
    Vector<N> toP = point-ray.p();
    if((ray.d()*toP)>0)
        return distanceSq(point,Line<N>(ray.p(), ray.d()));
    return toP.lengthSq();
}

template<unsigned int N>
double distanceSq(const Ray<N>& ray, const Vector<N>& point)
{
    return distanceSq(point,ray);
}

/**
* Computes Point to Ray Distance
* @param point
* @param ray
**/
template<unsigned int N>
double distance(const Vector<N>& point, const Ray<N>& ray)
{
    Vector<N> toP = point-ray.p();
    if((ray.d()*toP)>0)
        return distance(point,Line<N>(ray.p(), ray.d()));
    return toP.length();
}

template<unsigned int N>
double distance(const Ray<N>& ray, const Vector<N>& point)
{
    return distance(point,ray);
}

/**
* Computes Point to Segment Distance Squared
* @param point
* @param seg
**/
template<unsigned int N>
double distanceSq(const Vector<N>& point, const Segment<N>& seg)
{
    Vector<N> D = seg.d();
    Vector<N> toP = point-seg.p();
    double t = D*toP;
    if(t<=0)
    {
        return toP.lengthSq();
    }
    double DdD = D.length();
    if(t>=DdD)
    {
        Vector<N> toP1 = point-seg.P1();
        return toP1.lengthSq();
    }
    return distanceSq(point, Line<N>(seg.P0(), seg.d()));

}

template<unsigned int N>
double distanceSq(const Segment<N>& seg, const Vector<N>& point)
{
    return distanceSq(point,seg);
}

/**
* Computes Point to Segment Distance
* @param point
* @param seg
**/
template<unsigned int N>
double distance(const Vector<N>& point, const Segment<N>& seg)
{
    return sqrt(distanceSq(point,seg));
}

template<unsigned int N>
double distance(const Segment<N>& seg, const Vector<N>& point)
{
    return distance(point,seg);
}

} }

#endif
