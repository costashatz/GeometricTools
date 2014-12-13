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

#ifndef POLYGON_H
#define POLYGON_H

/**
* Includes
**/
#include <Primitives/2D/Polyline.h>

namespace GeometricTools {

namespace Primitives {

/**
* Polygon Class
* Polygon is a closed Polyline (closed means the first and the last points are the same)
* Assumes there is an extra point in the end identical to the first one
**/
class Polygon: public Polyline
{
public:
    /**
    * Default Constructor
    * Initialization
    **/
    Polygon():Polyline(){}

    /**
    * Get the area under the polygon
    * gives correct answer if polygon is simple/convex (=non-adjacent segments do not intersect)
    * @return double - the area
    **/
    virtual double area()
    {
        double sum = 0.0;
        for(int i=0;i<this->verts.size();i++)
        {
            unsigned int i_p = (i+1)%this->verts.size();
            sum += this->verts[i][0]*this->verts[i_p][1]-this->verts[i][1]*this->verts[i_p][0];
        }
        sum /= 2.0;
        return std::abs(sum);
    }

    /**
    * Check if vertices/points are clockwise ordered or not
    * @return bool - a boolean indicating if points are clockwise ordered
    **/
    bool clockwiseOrdered()
    {
        //clockwise (sum<0)
        double sum = 0.0;
        for(int i=0;i<this->verts.size();i++)
        {
            unsigned int i_p = (i+1)%this->verts.size();
            sum += this->verts[i][0]*this->verts[i_p][1]-this->verts[i][1]*this->verts[i_p][0];
        }
        if(sum < 0)
            return true;
        return false;
    }

    /**
    * Check if polygon is convex
    * @return bool - a boolean indicating whether the polygon is convex
    **/
    bool convex()
    {
        //convex (all cross products same sign)
        int plus=0,minus=0;
        for(int i=1;i<this->verts.size();i++)
        {
            int i_p = (i+1)%this->verts.size();
            if(((this->verts[i][0]-this->verts[i-1][0])*(this->verts[i_p][1]-this->verts[i][1])-(this->verts[i][1]-this->verts[i-1][1])*(this->verts[i_p][0]-this->verts[i][0]))<0)
                minus++;
            else
                plus++;
            if(plus>0&&minus>0)
                break;
        }
        if(plus>0&&minus>0)
            return false;
        return true;
    }
};

} }

#endif
