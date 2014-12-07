#ifndef QUADRATIC_CURVE_H
#define QUADRATIC_CURVE_H

/**
* Includes
**/
#include <Math/Matrix.h>


namespace LinearAlgebraTools {

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
    Matrix<2,2> A;
    Vector<2> B;
    double C;
public:
    QuadraticCurve() {}

    QuadraticCurve(const Matrix<2,2>& a, const Vector<2>& b, const double& c): A(a), B(b), C(c) {}

    Matrix<2,2> a() { return A; }

    Vector<2> b() { return B; }

    double c() { return C; }
};

} }

#endif
