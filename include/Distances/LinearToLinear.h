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

#ifndef LINEAR_TO_LINEAR_H
#define LINEAR_TO_LINEAR_H

/**
* Includes
**/
#include <Math/Vector.h>
#include <Primitives/LinearShapes.h>
#include <limits>

namespace GeometricTools {

using Math::Vector;
using Primitives::Line;
using Primitives::Ray;
using Primitives::Segment;

namespace Distances {

/**
* Computes Line to Line Distance Squared
* @param line1
* @param line2
**/
template<unsigned int N>
double distanceSq(const Line<N>& line1, const Line<N>& line2)
{
    Vector<N> w = line1.p()-line2.p();
    double a = line1.d()*line1.d();
    double b = line1.d()*line2.d();
    double c = line2.d()*line2.d();
    double d = line1.d()*w;
    double e = line2.d()*w;
    double D = a*c-b*b;
    double sC, tC;
    if (D<std::numeric_limits<double>::epsilon())
    {
        sC = 0.0;
        tC = (b>c)?(d/b):(e/c);
    }
    else
    {
        sC = (b*e-c*d)/D;
        tC = (a*e-b*d)/D;
    }

    Vector<N> dP = w+(sC*line1.d())-(tC*line2.d());
    return dP.lengthSq();
}

/**
* Computes Line to Line Distance
* @param line1
* @param line2
**/
template<unsigned int N>
double distance(const Line<N>& line1, const Line<N>& line2)
{
    return sqrt(distanceSq(line1,line2));
}

/**
* Computes Line to Ray Distance Squared
* @param line
* @param ray
**/
template<unsigned int N>
double distanceSq(const Line<N>& line, const Ray<N>& ray)
{
    Vector<N> u = line.d();
    Vector<N> v = ray.d();
    Vector<N> w = line.p()-ray.p();
    double a = u*u, b = u*v, c = v*v, d = u*w, e = v*w;
    double D = a*c-b*b;
    double sD = D, tD = D;
    double sN, tN, sc, tc;
    if(D<std::numeric_limits<double>::epsilon())
    {
        sN = 0.0;
        sD = 1.0;
        tN = e;
        tD = c;
    }
    else
    {
        sN = b*e-c*d;
        tN = a*e-b*d;
    }

    if(tN<0.0)
        tN = 0.0;

    if(std::abs(sN)<std::numeric_limits<double>::epsilon())
        sc = 0.0;
    else
        sc = sN/sD;

    if(std::abs(tN)<std::numeric_limits<double>::epsilon())
        tc = 0.0;
    else
        tc = tN/tD;

    Vector<N> dP = w+sc*u-tc*v;
    return dP.lengthSq();
}

template<unsigned int N>
double distanceSq(const Ray<N>& ray, const Line<N>& line)
{
    return distanceSq(line,ray);
}


/**
* Computes Line to Ray Distance
* @param line
* @param ray
**/
template<unsigned int N>
double distance(const Line<N>& line, const Ray<N>& ray)
{
    return sqrt(distanceSq(line,ray));
}

template<unsigned int N>
double distance(const Ray<N>& ray, const Line<N>& line)
{
    return distance(line,ray);
}

/**
* Computes Line to Segment Distance Squared
* @param line
* @param seg
**/
template<unsigned int N>
double distanceSq(const Line<N>& line, const Segment<N>& seg)
{
    Vector<N> u = line.d();
    Vector<N> v = seg.P0() - seg.P1();
    Vector<N> w = line.p() - seg.P0();
    double a = u*u, b = u*v, c = v*v, d = u*w, e = v*w;
    double D = a*c-b*b;
    double sD = D, tD = D;
    double sN, tN, sc, tc;
    if(D<std::numeric_limits<double>::epsilon())
    {
        sN = 0.0;
        sD = 1.0;
        tN = e;
        tD = c;
    }
    else
    {
        sN = b*e-c*d;
        tN = a*e-b*d;
    }

    if(tN<0.0)
    {
        tN = 0.0;
    }
    else if(tN>tD)
    {
        tN = tD;
    }

    if(std::abs(sN)<std::numeric_limits<double>::epsilon())
        sc = 0.0;
    else
        sc = sN/sD;

    if(std::abs(tN)<std::numeric_limits<double>::epsilon())
        tc = 0.0;
    else
        tc = tN/tD;

    Vector<N> dP = w+sc*u-tc*v;
    return dP.lengthSq();
}

template<unsigned int N>
double distanceSq(const Segment<N>& seg, const Line<N>& line)
{
    return distanceSq(line,seg);
}


/**
* Computes Line to Segment Distance
* @param line
* @param seg
**/
template<unsigned int N>
double distance(const Line<N>& line, const Segment<N>& seg)
{
    return sqrt(distanceSq(line,seg));
}

template<unsigned int N>
double distance(const Segment<N>& seg, const Line<N>& line)
{
    return distance(line,seg);
}

/**
* Computes Segment to Segment Distance Squared
* @param seg1
* @param seg2
**/
template<unsigned int N>
double distanceSq(const Segment<N>& seg1, const Segment<N>& seg2)
{
    Vector<N> u = seg1.P0() - seg1.P1();
    Vector<N> v = seg2.P0() - seg2.P1();
    Vector<N> w = seg1.P1() - seg2.P1();
    double a = u*u, b = u*v, c = v*v, d = u*w, e = v*w;
    double D = a*c-b*b;
    double sD = D, tD = D;
    double sN, tN, sc, tc;
    if(D<std::numeric_limits<double>::epsilon())
    {
        sN = 0.0;
        sD = 1.0;
        tN = e;
        tD = c;
    }
    else
    {
        sN = b*e-c*d;
        tN = a*e-b*d;
        if(sN<0.0)
        {
            sN = 0.0;
            tN = e;
            tD = c;
        }
        else if(sN>sD)
        {
            sN = sD;
            tN = e+b;
            tD = c;
        }
    }

    if(tN<0.0)
    {
        tN = 0.0;
        if(-d<0.0)
            sN = 0.0;
        else if(-d>a)
            sN = sD;
        else
        {
            sN = -d;
            sD = a;
        }
    }
    else if(tN>tD)
    {
        tN = tD;
        if((b-d)<0.0)
            sN = 0;
        else if((b-d)>a)
            sN = sD;
        else
        {
            sN = b-d;
            sD = a;
        }
    }

    if(std::abs(sN)<std::numeric_limits<double>::epsilon())
        sc = 0.0;
    else
        sc = sN/sD;

    if(std::abs(tN)<std::numeric_limits<double>::epsilon())
        tc = 0.0;
    else
        tc = tN/tD;

    Vector<N> dP = w+sc*u-tc*v;
    return dP.lengthSq();
}

/**
* Computes Segment to Segment Distance
* @param seg1
* @param seg2
**/
template<unsigned int N>
double distance(const Segment<N>& seg1, const Segment<N>& seg2)
{
    return sqrt(distanceSq(seg1, seg2));
}

} }


#endif
