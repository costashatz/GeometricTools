#ifndef LINEAR_TO_POLYLINE_H
#define LINEAR_TO_POLYLINE_H

/**
* Includes
**/
#include <Math/Vector.h>
#include <Primitives/LinearShapes.h>
#include <Distances/2D/PointToPolyline.h>
#include <Distances/2D/LinearToLinear.h>
#include <Primitives/2D/Polyline.h>

namespace LinearAlgebraTools {

using Math::Vector;
using Primitives::Line;
using Primitives::Ray;
using Primitives::Segment;
using Primitives::Polyline;

namespace Distances {


/**
* Computes Segment to Polyline Distance Squared
* @param segment
* @param polyline
**/
double DistanceSq(const Segment<2>& seg, const Polyline& polyline)
{
    double m = DistanceSq(seg, Segment<2>(polyline.vertices()[0], polyline.vertices()[1]));
    for(unsigned int i=1;i<polyline.vertices().size()-1;i++)
    {
        double tmp = DistanceSq(seg, Segment<2>(polyline.vertices()[i], polyline.vertices()[i+1]));
        if(tmp<m)
            m = tmp;
    }
    return m;
}

double DistanceSq(const Polyline& polyline, const Segment<2>& seg)
{
    return DistanceSq(seg, polyline);
}


/**
* Computes Segment to Polyline Distance
* @param segment
* @param polyline
**/
double Distance(const Segment<2>& seg, const Polyline& polyline)
{
    double m = Distance(seg, Segment<2>(polyline.vertices()[0], polyline.vertices()[1]));
    for(unsigned int i=1;i<polyline.vertices().size()-1;i++)
    {
        double tmp = Distance(seg, Segment<2>(polyline.vertices()[i], polyline.vertices()[i+1]));
        if(tmp<m)
            m = tmp;
    }
    return m;
}

double Distance(const Polyline& polyline, const Segment<2>& seg)
{
    return Distance(seg, polyline);
}

} }

#endif
