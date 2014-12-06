#ifndef POINT_TO_POLYLINE_H
#define POINT_TO_POLYLINE_H

/**
* Includes
**/
#include <Math/Vector.h>
#include <Primitives/LinearShapes.h>
#include <Distances/2D/PointToLinear.h>
#include <Primitives/2D/Polyline.h>

namespace LinearAlgebraTools {

using Math::Vector2;
using Primitives::Line2;
using Primitives::Ray2;
using Primitives::Segment2;
using Primitives::Polyline;

namespace Distances {

/**
* Computes Point to Polyline Distance Squared
* @param point
* @param line
**/
double DistanceSq(const Vector2& point, const Polyline& line)
{
    double m = DistanceSq(point, Segment2(line.vertices()[0], line.vertices()[1]));
    Vector2 prev = line.vertices()[1];
    double prevX = abs(prev[0]-point[0]);
    double prevY = abs(prev[1]-point[1]);
    double currX, currY;
    for(int i=1;i<line.vertices().size()-1;i++)
    {
        currX = abs(line.vertices()[i][0]-point[0]);
        currY = abs(line.vertices()[i][1]-point[1]);
        if((prevX <= m && currX <= m && (prevX*currX)>0) && (prevY <= m && currY <= m && (prevY*currY)>0))
        {
            m = DistanceSq(point, Segment2(line.vertices()[i], line.vertices()[i+1]));
        }
        prevX = currX;
        prevY = currY;
    }
    return m;
}

double DistanceSq(const Polyline& line, const Vector2& point)
{
    return DistanceSq(point,line);
}
/**
* Computes Point to Polyline Distance
* @param point
* @param line
**/
double Distance(const Vector2& point, const Polyline& line)
{
    double m = Distance(point, Segment2(line.vertices()[0], line.vertices()[1]));
    Vector2 prev = line.vertices()[1];
    double prevX = abs(prev[0]-point[0]);
    double prevY = abs(prev[1]-point[1]);
    double currX, currY;
    for(int i=1;i<line.vertices().size()-1;i++)
    {
        currX = abs(line.vertices()[i][0]-point[0]);
        currY = abs(line.vertices()[i][1]-point[1]);
        if((prevX <= m && currX <= m && (prevX*currX)>0) && (prevY <= m && currY <= m && (prevY*currY)>0))
        {
            m = Distance(point, Segment2(line.vertices()[i], line.vertices()[i+1]));
        }
        prevX = currX;
        prevY = currY;
    }
    return m;
}

double Distance(const Polyline& line, const Vector2& point)
{
    return Distance(point,line);
}

} }

#endif
