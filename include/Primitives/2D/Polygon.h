#ifndef POLYGON_H
#define POLYGON_H

#include <Primitives/2D/Polyline.h>

namespace LinearAlgebraTools {

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
    virtual double Area()
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
    bool ClockwiseOrdered()
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
    bool Convex()
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
