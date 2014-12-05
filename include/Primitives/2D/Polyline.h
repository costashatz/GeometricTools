#ifndef POLYLINE_H
#define POLYLINE_H

/**
* Includes
**/
#include <Math/Vector.h>
#include <vector>
#include <algorithm>
using std::vector;
using std::remove;

namespace LinearAlgebraTools {

using Math::Vector2;

namespace Primitives {

/**
* Polyline Class
* Polyline is a collection of arbitrary segments (or points)
* I use the form of collection of points..Edges are assumed to be P0->P1->P2->...->Pn
**/
class Polyline
{
protected:
    // Collection of points (called vertices)
    vector<Vector2> verts;
public:
    /**
    * Default Constructor
    * Initialization
    **/
    Polyline()
    {
        verts = vector<Vector2>();
    }

    /**
    * Add new point to the polyline
    * virtual method - can be overwritten by subclasses
    * @param point - point to be added
    **/
    virtual void AddPoint(const Vector2& point)
    {
        verts.push_back(point);
    }

    /**
    * Removes point from the polyline
    * virtual method - can be overwritten by subclasses
    * @param point - point to be removed
    **/
    virtual void RemovePoint(const Vector2& point)
    {
        remove(verts.begin(), verts.end(), point);
    }

    /**
    * Get Vertices/Points
    * @return vector<Vector2> - the collection of points/vertices
    **/
    vector<Vector2> vertices() {return verts;}
};

} }

#endif
