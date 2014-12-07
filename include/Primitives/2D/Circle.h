#ifndef CIRCLE_H
#define CIRCLE_H

/**
* Includes
**/
#include <Misc/Helper.h>
#include <Math/Matrix.h>
#include <Primitives/2D/QuadraticCurve.h>


namespace LinearAlgebraTools {

using Math::Vector;

namespace Primitives {

/**
* Circle Class
* Circle defined by center and radius
**/
class Circle : public QuadraticCurve
{
protected:
    Vector<2> Center;
    double Radius;
public:
    Circle(): QuadraticCurve() {}

    Circle(const Vector<2>& center, const double& radius): Center(center), Radius(radius)
    {
        A.identity();
        B = -2.0*center;
        C = center*center-radius*radius;
    }

    Vector<2> center() { return Center; }

    double radius() { return Radius; }

    double area() { return Helper::Pi*Radius*Radius; }
};

} }

#endif
