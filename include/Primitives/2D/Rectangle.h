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

#ifndef GEOMETRIC_TOOLS_PRIMITIVES_2D_RECTANGLE_H
#define GEOMETRIC_TOOLS_PRIMITIVES_2D_RECTANGLE_H

/**
* Includes
**/
#include <Math/Vector.h>
#include <Primitives/2D/Polygon.h>

namespace GeometricTools {

using Math::Vector2;

namespace Primitives {

/**
* Rectangle Class
**/
class Rectangle: public Polygon
{
public:
    /**
    * Default Constructor
    **/
    Rectangle(): Polygon() {}

    /**
    * Constructor
    * @param p0 - lower left point of rectangle
    * @param e0 - direction (with length) of horizontal edge
    * @param e1 - direction (with length) of vertical edge
    **/
    Rectangle(const Vector<2>& p0, const Vector<2>& e0, const Vector<2>& e1):Polygon()
    {
        Vector<2> E1 = e1;
        if((e0*e1)!=0)
            E1 = Vector<2>(-e0[0], e0[1]);
        Polygon::addPoint(p0);
        Polygon::addPoint(p0+e0);
        Polygon::addPoint(p0+e0+E1);
        Polygon::addPoint(p0+E1);
    }

    /**
    * Constructor
    * creates axis aligned rectangle
    * @param p - center of the rectangle
    * @param a - length of horizontal edge
    * @param b - length of vertical edge
    **/
    Rectangle(const Vector<2>& p, const double& a, const double& b)
    {
        Vector2 p0 = Vector2(p[0]-a/2.0, p[1]-b/2.0);
        Polygon::addPoint(p0);
        Vector2 p1 = Vector2(p[0]+a/2.0, p[1]-b/2.0);
        Polygon::addPoint(p1);
        Vector2 p2 = Vector2(p[0]+a/2.0, p[1]+b/2.0);
        Polygon::addPoint(p2);
        Vector2 p3 = Vector2(p[0]-a/2.0, p[1]+b/2.0);
        Polygon::addPoint(p3);
    }

    /**
    * Overwrite virtual method AddPoint so that it does nothing
    * We do not want other points to be added in a rectangle
    **/
    void addPoint(const Vector<2>& point)
    {
        return;
    }

    /**
    * Overwrite virtual method RemovePoint so that it does nothing
    * We do not want points to be removed in a rectangle
    **/
    void removePoint(const Vector<2>& point)
    {
        return;
    }
};

} }

#endif
