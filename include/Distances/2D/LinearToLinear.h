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
double DistanceSq(Line2 line1, Line2 line2)
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
double Distance(Line2 line1, Line2 line2)
{
    Vector2 d2P = Vector2(line2.d()[0], -line2.d()[1]);
    if(line1.d()*d2P!=0)
        return 0.0;
    Vector2 d1P = Vector2(line1.d()[0], -line1.d()[1]);
    double a = std::abs(d1P*(line2.p()-line1.p()));
    return a/line1.d().length();
}

} }


#endif
