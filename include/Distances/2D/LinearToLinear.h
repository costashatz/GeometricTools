#ifndef LINEAR_TO_LINEAR_H
#define LINEAR_TO_LINEAR_H

/**
* Includes
**/
#include <Math/Vector.h>
#include <Primitives/LinearShapes.h>

namespace LinearAlgebraTools {

using Math::Vector2;
using Primitives::Line2;
using Primitives::Ray2;
using Primitives::Segment2;

namespace Distances {

/**
* Computes Line to Line Distance Squared
* @param line1
* @param line2
**/
double DistanceSq(const Line2& line1, const Line2& line2)
{
    Vector2 d2P = Vector2(line2.d()[0], -line2.d()[1]);
    if(line1.d()*d2P!=0)
        return 0.0;
    Vector2 d1P = Vector2(line1.d()[0], -line1.d()[1]);
    double a = d1P*(line2.p()-line1.p());
    return (a*a)/line1.d().lengthSq();
}

/**
* Computes Line to Line Distance
* @param line1
* @param line2
**/
double Distance(const Line2& line1, const Line2& line2)
{
    Vector2 d2P = Vector2(line2.d()[0], -line2.d()[1]);
    if(line1.d()*d2P!=0)
        return 0.0;
    Vector2 d1P = Vector2(line1.d()[0], -line1.d()[1]);
    double a = std::abs(d1P*(line2.p()-line1.p()));
    return a/line1.d().length();
}

/**
* Computes Line to Ray Distance Squared
* @param line
* @param ray
**/
double DistanceSq(const Line2& line, const Ray2& ray)
{
    Vector2 d0 = Vector2(line.d()[1], -line.d()[0]);
    Vector2 d1 = ray.d();
    Vector2 D = ray.p()-line.p();
    if((d0*d1)*(d0*D)<0)
        return 0.0;
    double d = (d0*D);
    return d*d/d0.lengthSq();
}

double DistanceSq(const Ray2& ray, const Line2& line)
{
    return DistanceSq(line,ray);
}


/**
* Computes Line to Ray Distance
* @param line
* @param ray
**/
double Distance(const Line2& line, const Ray2& ray)
{
    Vector2 d0 = Vector2(line.d()[1], -line.d()[0]);
    Vector2 d1 = ray.d();
    Vector2 D = ray.p()-line.p();
    if((d0*d1)*(d0*D)<0)
        return 0.0;
    return std::abs(d0*D)/d0.length();
}

double Distance(const Ray2& ray, const Line2& line)
{
    return Distance(line,ray);
}

/**
* Computes Line to Segment Distance Squared
* @param line
* @param seg
**/
double DistanceSq(const Line2& line, const Segment2& seg)
{
    Vector2 d0 = Vector2(line.d()[1], -line.d()[0]);
    Vector2 d1 = seg.d();
    double T1 = 1.0;
    Vector2 D = seg.p()-line.p();
    if(((d0*D)*(d0*(D+T1*d1)))<0)
        return 0.0;
    double d = std::min(std::abs(d0*D), std::abs(d0*(D+T1*d1)));
    return d*d/d0.lengthSq();
}

double DistanceSq(const Segment2& seg, const Line2& line)
{
    return DistanceSq(line,seg);
}


/**
* Computes Line to Segment Distance
* @param line
* @param seg
**/
double Distance(const Line2& line, const Segment2& seg)
{
    Vector2 d0 = Vector2(line.d()[1], -line.d()[0]);
    Vector2 d1 = seg.d();
    double T1 = 1.0;
    Vector2 D = seg.p()-line.p();
    if(((d0*D)*(d0*(D+T1*d1)))<0)
        return 0.0;
    return std::min(std::abs(d0*D), std::abs(d0*(D+T1*d1)))/d0.length();
}

double Distance(const Segment2& seg, const Line2& line)
{
    return Distance(line,seg);
}

} }


#endif
