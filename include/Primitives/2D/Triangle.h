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

#ifndef GEOMETRIC_TOOLS_PRIMITIVES_2D_TRIANGLE_H
#define GEOMETRIC_TOOLS_PRIMITIVES_2D_TRIANGLE_H

/**
* Includes
**/
#include <Primitives/2D/Polygon.h>

namespace GeometricTools {

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
    Triangle(const Vector<2>& p0, const Vector<2>& p1, const Vector<2>& p2):Polygon()
    {
        Polygon::addPoint(p0);
        Polygon::addPoint(p1);
        Polygon::addPoint(p2);
    }

    /**
    * Overwrite virtual method AddPoint so that it does nothing
    * We do not want other points to be added in a triangle
    **/
    void addPoint(const Vector<2>& point)
    {
        return;
    }

    /**
    * Overwrite virtual method RemovePoint so that it does nothing
    * We do not want points to be removed in a triangle
    **/
    void removePoint(const Vector<2>& point)
    {
        return;
    }
};

} }

#endif
