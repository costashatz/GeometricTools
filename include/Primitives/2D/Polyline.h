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

namespace GeometricTools {

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
