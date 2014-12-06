#ifndef RECTANGLE_H
#define RECTANGLE_H

/**
* Includes
**/
#include <Math/Vector.h>
#include <Primitives/2D/Polygon.h>

namespace LinearAlgebraTools {

using Math::Vector2;

namespace Primitives {

/**
* Rectangle Class
**/
template<class T>
class Rectangle: public Polygon
{
public:
    /**
    * Constructor
    * @param p0 - lower left point of rectangle
    * @param e0 - direction (with length) of horizontal edge
    * @param e1 - direction (with length) of vertical edge
    **/
    Rectangle(const Vector2& p0, const Vector2& e0, const Vector2& e1):Polygon()
    {
        if((e0*e1)!=0)
            e1 = Vector2(-e0[0], e0[1]);
        Polygon::AddPoint(p0);
        Polygon::AddPoint(p0+e0);
        Polygon::AddPoint(p0+e0+e1);
        Polygon::AddPoint(p0+e1);
    }

    /**
    * Constructor
    * creates axis aligned rectangle
    * @param p - center of the rectangle
    * @param a - length of horizontal edge
    * @param b - length of vertical edge
    **/
    Rectangle(const Vector2& p, const double& a, const double& b)
    {
        Vector2 p0 = Vector2(p[0]-a/2.0, p[1]-b/2.0);
        Polygon::AddPoint(p0);
        Vector2 p1 = Vector2(p[0]+a/2.0, p[1]-b/2.0);
        Polygon::AddPoint(p1);
        Vector2 p2 = Vector2(p[0]+a/2.0, p[1]+b/2.0);
        Polygon::AddPoint(p2);
        Vector2 p3 = Vector2(p[0]-a/2.0, p[1]+b/2.0);
        Polygon::AddPoint(p3);
    }

    /**
    * Overwrite virtual method AddPoint so that it does nothing
    * We do not want other points to be added in a triangle
    **/
    void AddPoint(const Vector2& point)
    {
        return;
    }

    /**
    * Overwrite virtual method RemovePoint so that it does nothing
    * We do not want points to be removed in a triangle
    **/
    void RemovePoint(const Vector2& point)
    {
        return;
    }
};

} }

#endif
