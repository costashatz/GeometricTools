#ifndef TRIANGLE_H
#define TRIANGLE_H

#include <Primitives/2D/Polygon.h>

namespace LinearAlgebraTools {

using Math::Vector2;

namespace Primitives {

/**
* Triangle Class
**/
class Triangle: public Polygon
{
public:
    /**
    * Constructor
    * @param p0 - the first point of the triangle
    * @param p1 - the second point of the triangle
    * @param p2 - the third point of the triangle
    **/
    Triangle(const Vector2& p0, const Vector2& p1, const Vector2& p2):Polygon()
    {
        Polygon::AddPoint(p0);
        Polygon::AddPoint(p1);
        Polygon::AddPoint(p2);
    }

    /**
    * Overwrite virtual method AddPoint so that it does nothing
    * We do not want other points to be added in a triangle
    **/
    void AddPoint(const Vector2& point)
    {
        return;
    }
};

} }

#endif
