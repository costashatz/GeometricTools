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

using Math::Vector;

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
    vector<Vector<2> > verts;
public:
    /**
    * Default Constructor
    * Initialization
    **/
    Polyline()
    {
        verts = vector<Vector<2> >();
    }

    /**
    * Add new point to the polyline
    * virtual method - can be overwritten by subclasses
    * @param point - point to be added
    **/
    virtual void addPoint(const Vector<2>& point)
    {
        verts.push_back(point);
    }

    /**
    * Removes point from the polyline
    * virtual method - can be overwritten by subclasses
    * @param point - point to be removed
    **/
    virtual void removePoint(const Vector<2>& point)
    {
        remove(verts.begin(), verts.end(), point);
    }

    /**
    * Get Vertices/Points
    * @return vector<Vector<2> > - the collection of points/vertices
    **/
    vector<Vector<2> > vertices() const {return verts;}
};

} }

#endif
