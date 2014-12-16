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

#ifndef GEOMETRIC_TOOLS_PRIMITIVES_TOOLS_BOUNDING_BOX_H
#define GEOMETRIC_TOOLS_PRIMITIVES_TOOLS_BOUNDING_BOX_H

/**
* Includes
**/
#include <Primitives/2D/Rectangle.h>
#include <Primitives/2D/Polygon.h>
#include <Math/Vector.h>
#include <vector>

using std::vector;

namespace GeometricTools {

using Math::Vector;

namespace Primitives {

inline Rectangle boundingBox(const vector<Vector<2> >& points)
{
    //TODO: assert points size > 1
    double minX, maxX, minY, maxY;
    minX = maxX = points[0][0];
    minY = maxY = points[0][1];

    for(int i=1;i<points.size();i++)
    {
        if(points[i][0]<minX)
            minX = points[i][0];
        else if(points[i][0]>maxX)
            maxX = points[i][0];
        if(points[i][1]<minY)
            minY = points[i][1];
        else if(points[i][1]>maxY)
            maxY = points[i][1];
    }

    double a = maxX-minX;
    double b = maxY-minY;
    return Rectangle(Vector<2>(minX+a/2, minY+b/2), a, b);
}

inline Rectangle boundingBox(const vector<Polygon>& polygons)
{
    //TODO: assert polygons size >= 1
    vector<Vector<2> > points = polygons[0].vertices();
    for(int i=1;i<polygons.size();i++)
    {
        for(int j=0;j<polygons[i].vertices().size();j++)
            points.push_back(polygons[i].vertices()[j]);
    }
    return boundingBox(points);
}

} }

#endif
