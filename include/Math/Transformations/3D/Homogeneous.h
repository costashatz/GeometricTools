#ifndef HOMOGENEOUS_3D_H
#define HOMOGENEOUS_3D_H

/**
* Includes
**/
#include <Math/Matrix.h>
#include <Misc/Helper.h>
#include <cmath>

namespace LinearAlgebraTools { namespace Math {

namespace Transformations3D {

/**
* Get 3D translation matrix
* @param dx - x translation
* @param dy - y translation
* @param dz - z translation
* @return Matrix<4,4> - translation matrix
**/
Matrix<4,4> translation(const double& dx, const double& dy, const double& dz)
{
    return Matrix<4,4>(1.0, 0.0, 0.0, dx,
                       0.0, 1.0, 0.0, dy,
                       0.0, 0.0, 1.0, dz,
                       0.0, 0.0, 0.0, 1.0);
}

/**
* Get 3D scaling matrix
* @param sx - x scaling factor
* @param sy - y scaling factor
* @param sz - z scaling factor
* @return Matrix<4,4> - scaling matrix
**/
Matrix<4,4> scaling(const double& sx, const double& sy, const double& sz)
{
    return Matrix<4,4>(sx , 0.0, 0.0, 0.0,
                       0.0, sy , 0.0, 0.0,
                       0.0, 0.0, sz , 0.0,
                       0.0, 0.0, 0.0, 1.0);
}

/**
* Get 3D uniform scaling matrix
* @param s - scaling factor
* @return Matrix<4,4> - uniform scaling matrix
**/
Matrix<4,4> uniformScaling(const double& s)
{
    return Matrix<4,4>( s , 0.0, 0.0, 0.0,
                       0.0,  s , 0.0, 0.0,
                       0.0, 0.0,  s , 0.0,
                       0.0, 0.0, 0.0, 1.0);
}

/**
* Get 3D shear matrix in x
* @param a - y shearing factor
* @param b - z shearing factor
* @return Matrix<4,4> - shear matrix in x
**/
Matrix<4,4> shearX(const double& a, const double& b)
{
    return Matrix<4,4>(1.0, 0.0, 0.0, 0.0,
                        a , 1.0, 0.0, 0.0,
                        b , 0.0, 1.0, 0.0,
                       0.0, 0.0, 0.0, 1.0);
}

/**
* Get 3D shear matrix in y
* @param a - x shearing factor
* @param b - z shearing factor
* @return Matrix<4,4> - shear matrix in y
**/
Matrix<4,4> shearY(const double& a, const double& b)
{
    return Matrix<4,4>(1.0,  a , 0.0, 0.0,
                       0.0, 1.0, 0.0, 0.0,
                       0.0,  b , 1.0, 0.0,
                       0.0, 0.0, 0.0, 1.0);
}

/**
* Get 3D shear matrix in z
* @param a - x shearing factor
* @param b - y shearing factor
* @return Matrix<4,4> - shear matrix in z
**/
Matrix<4,4> shearZ(const double& a, const double& b)
{
    return Matrix<4,4>(1.0, 0.0,  a , 0.0,
                       0.0, 1.0,  b , 0.0,
                       0.0, 0.0, 1.0, 0.0,
                       0.0, 0.0, 0.0, 1.0);
}

/**
* Get 3D rotation matrix around x
* @param a - angle in radians
* @return Matrix<4,4> - rotation matrix around x
**/
Matrix<4,4> rotationX(const double& a)
{
    return Matrix<4,4>(1.0,   0.0  ,  0.0  , 0.0,
                       0.0,  cos(a), sin(a), 0.0,
                       0.0, -sin(a), cos(a), 0.0,
                       0.0,   0.0  ,  0.0  , 1.0);
}

/**
* Get 3D rotation matrix around x
* @param a - angle in degrees
* @return Matrix<4,4> - rotation matrix around x
**/
Matrix<4,4> rotationXDegrees(const double& a)
{
    return rotationX(Helper::degreesToRadians(a));
}

/**
* Get 3D rotation matrix around y
* @param a - angle in radians
* @return Matrix<4,4> - rotation matrix around y
**/
Matrix<4,4> rotationY(const double& a)
{
    return Matrix<4,4>(cos(a),  0.0, sin(a), 0.0,
                       0.0,     1.0,  0.0  , 0.0,
                       -sin(a), 0.0, cos(a), 0.0,
                       0.0,   0.0  ,  0.0  , 1.0);
}

/**
* Get 3D rotation matrix around y
* @param a - angle in degrees
* @return Matrix<4,4> - rotation matrix around y
**/
Matrix<4,4> rotationYDegrees(const double& a)
{
    return rotationY(Helper::degreesToRadians(a));
}

/**
* Get 3D rotation matrix around z
* @param a - angle in radians
* @return Matrix<4,4> - rotation matrix around z
**/
Matrix<4,4> rotationZ(const double& a)
{
    return Matrix<4,4>( cos(a), sin(a), 0.0, 0.0,
                       -sin(a), cos(a), 0.0, 0.0,
                         0.0  ,  0.0  , 1.0, 0.0,
                         0.0  ,  0.0  , 0.0, 1.0);
}

/**
* Get 3D rotation matrix around z
* @param a - angle in degrees
* @return Matrix<4,4> - rotation matrix around z
**/
Matrix<4,4> rotationZDegrees(const double& a)
{
    return rotationZ(Helper::degreesToRadians(a));
}

/**
* Get 3D rotation matrix around axis
* @param axis - Vector3 axis to rotate around
* @param a - angle in radians
* @return Matrix<4,4> - rotation matrix around z
**/
Matrix<4,4> rotationAxis(const Vector3& axis, const double& a)
{
    double ux_2 = axis[0]*axis[0], uxuy = axis[0]*axis[1], uxuz = axis[0]*axis[2], uyuz = axis[1]*axis[2];
    double c = cos(a), s = sin(a), c_1 = 1.0-c, uy_2 = axis[1]*axis[1], uz_2 = axis[2]*axis[2];
    double ux = axis[0], uy = axis[1], uz = axis[2];
    return Matrix<4,4>(c+ux_2*c_1, uxuy*c_1-uz*s, uxuz*c_1+uy*s, 0.0,
                       uxuy*c_1+uz*s, c+uy_2*c_1, uyuz*c_1-ux*s, 0.0,
                       uxuz*c_1-uy*s, uyuz*c_1+ux*s, c+uz_2*c_1, 0.0,
                            0.0     ,      0.0     ,     0.0   , 1.0);
}

/**
* Get 3D rotation matrix around axis
* @param axis - Vector3 axis to rotate around
* @param a - angle in degrees
* @return Matrix<4,4> - rotation matrix around z
**/
Matrix<4,4> rotationAxisDegrees(const Vector3& axis, const double& a)
{
    return rotationAxis(axis, Helper::degreesToRadians(a));
}

} } }

#endif
