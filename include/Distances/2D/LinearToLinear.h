#ifndef LINEAR_TO_LINEAR_H
#define LINEAR_TO_LINEAR_H

/**
* Includes
**/
#include <Math/Vector.h>
#include <Primitives/LinearShapes.h>
#include <limits>

namespace LinearAlgebraTools {

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
double DistanceSq(const Line<2>& line1, const Line<2>& line2)
{
    Vector<2> d2P = Vector<2>(line2.d()[0], -line2.d()[1]);
    if(line1.d()*d2P!=0)
        return 0.0;
    Vector<2> d1P = Vector<2>(line1.d()[0], -line1.d()[1]);
    double a = d1P*(line2.p()-line1.p());
    return (a*a)/line1.d().lengthSq();
}

/**
* Computes Line to Line Distance
* @param line1
* @param line2
**/
double Distance(const Line<2>& line1, const Line<2>& line2)
{
    Vector<2> d2P = Vector<2>(line2.d()[0], -line2.d()[1]);
    if(line1.d()*d2P!=0)
        return 0.0;
    Vector<2> d1P = Vector<2>(line1.d()[0], -line1.d()[1]);
    double a = std::abs(d1P*(line2.p()-line1.p()));
    return a/line1.d().length();
}

/**
* Computes Line to Ray Distance Squared
* @param line
* @param ray
**/
double DistanceSq(const Line<2>& line, const Ray<2>& ray)
{
    Vector<2> d0 = Vector<2>(line.d()[1], -line.d()[0]);
    Vector<2> d1 = ray.d();
    Vector<2> D = ray.p()-line.p();
    if((d0*d1)*(d0*D)<0)
        return 0.0;
    double d = (d0*D);
    return d*d/d0.lengthSq();
}

double DistanceSq(const Ray<2>& ray, const Line<2>& line)
{
    return DistanceSq(line,ray);
}


/**
* Computes Line to Ray Distance
* @param line
* @param ray
**/
double Distance(const Line<2>& line, const Ray<2>& ray)
{
    Vector<2> d0 = Vector<2>(line.d()[1], -line.d()[0]);
    Vector<2> d1 = ray.d();
    Vector<2> D = ray.p()-line.p();
    if((d0*d1)*(d0*D)<0)
        return 0.0;
    return std::abs(d0*D)/d0.length();
}

double Distance(const Ray<2>& ray, const Line<2>& line)
{
    return Distance(line,ray);
}

/**
* Computes Line to Segment Distance Squared
* @param line
* @param seg
**/
double DistanceSq(const Line<2>& line, const Segment<2>& seg)
{
    Vector<2> d0 = Vector<2>(line.d()[1], -line.d()[0]);
    Vector<2> d1 = seg.d();
    double T1 = 1.0;
    Vector<2> D = seg.p()-line.p();
    if(((d0*D)*(d0*(D+T1*d1)))<0)
        return 0.0;
    double d = std::min(std::abs(d0*D), std::abs(d0*(D+T1*d1)));
    return d*d/d0.lengthSq();
}

double DistanceSq(const Segment<2>& seg, const Line<2>& line)
{
    return DistanceSq(line,seg);
}


/**
* Computes Line to Segment Distance
* @param line
* @param seg
**/
double Distance(const Line<2>& line, const Segment<2>& seg)
{
    Vector<2> d0 = Vector<2>(line.d()[1], -line.d()[0]);
    Vector<2> d1 = seg.d();
    double T1 = 1.0;
    Vector<2> D = seg.p()-line.p();
    if(((d0*D)*(d0*(D+T1*d1)))<0)
        return 0.0;
    return std::min(std::abs(d0*D), std::abs(d0*(D+T1*d1)))/d0.length();
}

double Distance(const Segment<2>& seg, const Line<2>& line)
{
    return Distance(line,seg);
}

/**
* Computes Segment to Segment Distance Squared
* @param seg1
* @param seg2
**/
double DistanceSq(const Segment<2>& seg1, const Segment<2>& seg2)
{
    Vector<2> u = seg1.P0() - seg1.P1();
    Vector<2> v = seg2.P0() - seg2.P1();
    Vector<2> w = seg1.P1() - seg2.P1();
    double a = u*u, b = u*v, c = v*v, d = u*w, e = v*w;
    double D = a*c-b*b;
    double sD = D, tD = D;
    double sN, tN, sc, tc;
    if(std::abs(D)<std::numeric_limits<double>::epsilon())
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

    Vector<2> dP = w+sc*u-tc*v;
    return dP.lengthSq();
}

/**
* Computes Segment to Segment Distance
* @param seg1
* @param seg2
**/
double Distance(const Segment<2>& seg1, const Segment<2>& seg2)
{
    return sqrt(DistanceSq(seg1, seg2));
}

} }


#endif
