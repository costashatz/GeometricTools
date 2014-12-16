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

#ifndef GEOMETRIC_TOOLS_PRIMITIVES_LINEAR_SHAPES_H
#define GEOMETRIC_TOOLS_PRIMITIVES_LINEAR_SHAPES_H

/**
* Includes
**/
#include <Math/Vector.h>

namespace GeometricTools {

using Math::Vector;

namespace Primitives {

/**
* Base Class for Linear Shapes
* I have used the parametric form of linear shape: X(t) = P+t*d, P=point, d=direction vector, t=parameter
**/
template<unsigned int N>
class LinearShape
{
protected:
    // LinearShape data
    Vector<N> P;
    Vector<N> D;
public:
    /**
    * Default Constructor
    * Initializes both p and d to zero vectors
    **/
    LinearShape() {}

    /**
    * Constructor - Parametric Form
    * @param p - point P
    * @param d - direction vector d
    **/
    LinearShape(const Vector<N>& p, const Vector<N>& d):P(p),D(d){}

    /**
    * Copy Constructor
    * @param other - LinearShape to copy from
    **/
    LinearShape(const LinearShape& other)
    {
        P = other.P;
        D = other.D;
    }

    /**
    * Get Point P - Parametric Form
    * @return Vector - the Point p
    **/
    Vector<N> p()const {return P;}

    /**
    * Get Direction Vector d - Parametric Form
    * @return Vector - the direction vector
    **/
    Vector<N> d()const {return D;}
};

/**
* Line Class
* A line is a linear shape with no restrictions on "t" value of the parametric form
**/
template<unsigned int N>
class Line: public LinearShape<N>
{
public:
    /**
    * Default Constructor
    **/
    Line():LinearShape<N>(){}

    /**
    * Constructor - Parametric Form
    * @param P - point P
    * @param D - direction vector d
    **/
    Line(const Vector<N>& P, const Vector<N>& D):LinearShape<N>(P,D)
    {
        this->D.normalize();
    }
};


/**
* Ray Class
* A ray is a linear shape with t>=0 (parametric form)
**/
template<unsigned int N>
class Ray: public LinearShape<N>
{
public:
    /**
    * Default Constructor
    **/
    Ray():LinearShape<N>(){}

    /**
    * Constructor - A ray is defined by a starting point P and a direction D
    * @param P - point P
    * @param D - direction vector d
    **/
    Ray(const Vector<N>& P, const Vector<N>& D):LinearShape<N>(P,D)
    {
        this->D.normalize();
    }
};


/**
* Segment Class
* A segment is a linear shape with tE[0,1] (parametric form)
**/
template<unsigned int N>
class Segment: public LinearShape<N>
{
public:
    /**
    * Default Constructor
    **/
    Segment():LinearShape<N>(){}

    /**
    * Constructor - A segment is defined by two points (P0-starting point, P1-ending point)
    * @param P0 - starting point
    * @param P1 - ending point
    **/
    Segment(const Vector<N>& P0, const Vector<N>& P1):LinearShape<N>(P0, P1-P0){}


    /**
    * Get Starting Point P0 - Is exactly the same as P() - I've included it for clarity/completeness
    * @return Vector - the Starting Point P0
    **/
    Vector<N> P0()const {return this->P;}


    /**
    * Get Ending Point P1
    * @return Vector - the Ending Point P1
    **/
    Vector<N> P1()const {return this->P+this->D;}

    /**
    * Get Length of Segment
    * @return double - length of segment
    **/
    double length()const { return this->D.length();}

    /**
    * Get Length Squared of Segment
    * @return double - length squared of segment
    **/
    double lengthSq()const { return this->D.lengthSq();}
};

/**
* Typedefs for frequently used types
**/
typedef LinearShape<2> LinearShape2;
typedef LinearShape<3> LinearShape3;
typedef Line<2> Line2;
typedef Line<3> Line3;
typedef Ray<2> Ray2;
typedef Ray<3> Ray3;
typedef Segment<2> Segment2;
typedef Segment<3> Segment3;

} }

#endif
