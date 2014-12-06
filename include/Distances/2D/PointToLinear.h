#ifndef POINT_TO_LINEAR_H
#define POINT_TO_LINEAR_H

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
* Computes Point to Line Distance Squared
* @param point
* @param line
**/
double DistanceSq(Vector2 point, Line2 line)
{
    Vector2 y_p = point-line.p();
    double k = (line.d()*y_p);
    return (y_p.lengthSq()-k*k);
}

double DistanceSq(Line2 line, Vector2 point)
{
    return DistanceSq(point,line);
}

/**
* Computes Point to Line Distance
* @param point
* @param line
**/
double Distance(Vector2 point, Line2 line)
{
    return sqrt(DistanceSq(point, line));
}

double Distance(Line2 line, Vector2 point)
{
    return Distance(point,line);
}

/**
* Computes Point to Ray Distance Squared
* @param point
* @param ray
**/
double DistanceSq(Vector2 point, Ray2 ray)
{
    Vector2 toP = point-ray.p();
    if((ray.d()*toP)>0)
        return DistanceSq(point,Line2(ray.p(), ray.d()));
    return toP.lengthSq();
}

double DistanceSq(Ray2 ray, Vector2 point)
{
    return DistanceSq(point,ray);
}

/**
* Computes Point to Ray Distance
* @param point
* @param ray
**/
double Distance(Vector2 point, Ray2 ray)
{
    Vector2 toP = point-ray.p();
    if((ray.d()*toP)>0)
        return Distance(point,Line2(ray.p(), ray.d()));
    return toP.length();
}

double Distance(Ray2 ray, Vector2 point)
{
    return Distance(point,ray);
}

/**
* Computes Point to Segment Distance Squared
* @param point
* @param seg
**/
double DistanceSq(Vector2 point, Segment2 seg)
{
    Vector2 D = seg.d();
    Vector2 toP = point-seg.p();
    double t = D*toP;
    if(t<=0)
    {
        return toP.lengthSq();
    }
    double DdD = D.lengthSq();
    if(t>=DdD)
    {
        Vector2 toP1 = point-seg.P1();
        return toP1.lengthSq();
    }
    return DistanceSq(point, Line2(seg.P0(), seg.d()));
}

double DistanceSq(Segment2 seg, Vector2 point)
{
    return DistanceSq(point,seg);
}

/**
* Computes Point to Segment Distance
* @param point
* @param seg
**/
double Distance(Vector2 point, Segment2 seg)
{
    Vector2 D = seg.d();
    Vector2 toP = point-seg.p();
    double t = D*toP;
    if(t<=0)
    {
        return toP.length();
    }
    double DdD = D.lengthSq();
    if(t>=DdD)
    {
        Vector2 toP1 = point-seg.P1();
        return toP1.length();
    }
    return Distance(point, Line2(seg.P0(), seg.d()));
}

double Distance(Segment2 seg, Vector2 point)
{
    return Distance(point,seg);
}

} }

#endif
