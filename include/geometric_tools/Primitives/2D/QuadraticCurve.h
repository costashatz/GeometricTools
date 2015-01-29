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

#ifndef GEOMETRIC_TOOLS_PRIMITIVES_2D_QUADRATIC_CURVE_H
#define GEOMETRIC_TOOLS_PRIMITIVES_2D_QUADRATIC_CURVE_H

/**
* Includes
**/
#include <geometric_tools/Math/Matrix.h>


namespace GeometricTools {

using Math::Matrix;
using Math::Vector;

namespace Primitives {

/**
* QuadraticCurve Class
* Quadratic curve defined by xTAx+bTx+c=0
**/
class QuadraticCurve
{
protected:
    Matrix<2,2> a_;
    Vector<2> b_;
    double c_;
public:
    QuadraticCurve() {}

    QuadraticCurve(const Matrix<2,2>& a, const Vector<2>& b, const double& c): a_(a), b_(b), c_(c) {}

    Matrix<2,2> a() { return a_; }

    Vector<2> b() { return b_; }

    double c() { return c_; }
};

} }

#endif
