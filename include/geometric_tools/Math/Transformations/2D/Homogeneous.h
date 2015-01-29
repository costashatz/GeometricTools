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

#ifndef GEOMETRIC_TOOLS_MATH_TRANSFORMATIONS_2D_HOMOGENEOUS_H
#define GEOMETRIC_TOOLS_MATH_TRANSFORMATIONS_2D_HOMOGENEOUS_H

/**
* Includes
**/
#include <geometric_tools/Math/Matrix.h>
#include <geometric_tools/Misc/Helper.h>
#include <cmath>

namespace GeometricTools { namespace Math {

namespace Transformations2D {

/**
* Get 2D translation matrix
* @param dx - x translation
* @param dy - y translation
* @return Matrix<3,3> - translation matrix
**/
inline Matrix<3,3> translation(const double& dx, const double& dy)
{
    return Matrix<3,3>(1.0, 0.0, dx,
                       0.0, 1.0, dy,
                       0.0, 0.0, 1.0);
}

/**
* Get 2D scaling matrix
* @param sx - x scaling factor
* @param sy - y scaling factor
* @return Matrix<3,3> - scaling matrix
**/
inline Matrix<3,3> scaling(const double& sx, const double& sy)
{
    return Matrix<3,3>(sx , 0.0, 0.0,
                       0.0, sy , 0.0,
                       0.0, 0.0, 1.0);
}

/**
* Get 2D uniform scaling matrix
* @param s - scaling factor
* @return Matrix<3,3> - uniform scaling matrix
**/
inline Matrix<3,3> uniformScaling(const double& s)
{
    return Matrix<3,3>( s , 0.0, 0.0,
                       0.0,  s , 0.0,
                       0.0, 0.0, 1.0);
}

/**
* Get 2D shear matrix
* @param sx - x shearing factor
* @param sy - y shearing factor
* @return Matrix<3,3> - shear matrix
**/
inline Matrix<3,3> shear(const double& sx, const double& sy)
{
    return Matrix<3,3>(1.0, sx , 0.0,
                       sy , 0.0, 0.0,
                       0.0, 0.0, 1.0);
}

/**
* Get 2D rotation matrix
* @param a - angle in radians
* @return Matrix<3,3> - rotation matrix
**/
inline Matrix<3,3> rotation(const double& a)
{
    return Matrix<3,3>( cos(a), -sin(a), 0.0,
                        sin(a),  cos(a), 0.0,
                        0.0,     0.0,    1.0);
}

/**
* Get 2D rotation matrix
* @param a - angle in degrees
* @return Matrix<3,3> - rotation matrix
**/
inline Matrix<3,3> rotationDegrees(const double& a)
{
    return rotation(Helper::degreesToRadians(a));
}

} } }

#endif
