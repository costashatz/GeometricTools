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

#ifndef GEOMETRIC_TOOLS_DISTANCES_2D_POINT_TO_POLYLINE_H
#define GEOMETRIC_TOOLS_DISTANCES_2D_POINT_TO_POLYLINE_H

/**
* Includes
**/
#include <Math/Vector.h>
#include <Primitives/LinearShapes.h>
#include <Distances/PointToLinear.h>
#include <Primitives/Polyline.h>

namespace GeometricTools {

using Math::Vector;
using Primitives::Line;
using Primitives::Ray;
using Primitives::Segment;
using Primitives::Polyline;

namespace Distances {

/**
* Computes Point to Polyline Distance Squared
* @param point
* @param line
**/
inline double distanceSq(const Vector<2>& point, const Polyline<2>& line)
{
    double m = distanceSq(point, Segment<2>(line.vertices()[0], line.vertices()[1]));
    Vector<2> prev = line.vertices()[1];
    double prevX = abs(prev[0]-point[0]);
    double prevY = abs(prev[1]-point[1]);
    double currX, currY;
    for(int i=1;i<line.vertices().size()-1;i++)
    {
        currX = abs(line.vertices()[i][0]-point[0]);
        currY = abs(line.vertices()[i][1]-point[1]);
        if((prevX <= m && currX <= m && (prevX*currX)>0) && (prevY <= m && currY <= m && (prevY*currY)>0))
        {
            m = distanceSq(point, Segment<2>(line.vertices()[i], line.vertices()[i+1]));
        }
        prevX = currX;
        prevY = currY;
    }
    return m;
}

inline double distanceSq(const Polyline<2>& line, const Vector<2>& point)
{
    return distanceSq(point,line);
}
/**
* Computes Point to Polyline Distance
* @param point
* @param line
**/
inline double distance(const Vector<2>& point, const Polyline<2>& line)
{
    double m = distance(point, Segment<2>(line.vertices()[0], line.vertices()[1]));
    Vector<2> prev = line.vertices()[1];
    double prevX = abs(prev[0]-point[0]);
    double prevY = abs(prev[1]-point[1]);
    double currX, currY;
    for(int i=1;i<line.vertices().size()-1;i++)
    {
        currX = abs(line.vertices()[i][0]-point[0]);
        currY = abs(line.vertices()[i][1]-point[1]);
        if((prevX <= m && currX <= m && (prevX*currX)>0) && (prevY <= m && currY <= m && (prevY*currY)>0))
        {
            m = distance(point, Segment<2>(line.vertices()[i], line.vertices()[i+1]));
        }
        prevX = currX;
        prevY = currY;
    }
    return m;
}

inline double distance(const Polyline<2>& line, const Vector<2>& point)
{
    return distance(point,line);
}

} }

#endif
