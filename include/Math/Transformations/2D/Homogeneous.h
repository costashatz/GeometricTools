#ifndef HOMOGENEOUS_2D_H
#define HOMOGENEOUS_2D_H

/**
* Includes
**/
#include <Math/Matrix.h>
#include <Misc/Helper.h>
#include <cmath>

namespace LinearAlgebraTools { namespace Math {

namespace Transformations2D {

/**
* Get 2D translation matrix
* @param dx - x translation
* @param dy - y translation
* @return Matrix<3,3> - translation matrix
**/
Matrix<3,3> translation(const double& dx, const double& dy)
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
Matrix<3,3> scaling(const double& sx, const double& sy)
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
Matrix<3,3> uniformScaling(const double& s)
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
Matrix<3,3> shear(const double& sx, const double& sy)
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
Matrix<3,3> rotation(const double& a)
{
    return Matrix<3,3>( cos(a), sin(a), 0.0,
                       -sin(a), cos(a), 0.0,
                        0.0,    0.0,    1.0);
}

/**
* Get 2D rotation matrix
* @param a - angle in degrees
* @return Matrix<3,3> - rotation matrix
**/
Matrix<3,3> rotationDegrees(const double& a)
{
    return rotation(Helper::degreesToRadians(a));
}

} } }

#endif
