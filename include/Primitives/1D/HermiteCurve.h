#ifndef HERMITE_CURVE_H
#define HERMITE_CURVE_H

/**
* Includes
**/
#include <Math/Vector.h>
#include <Math/Matrix.h>
#include <Math/LinearSystems/SolveLU.h>
#include <Primitives/1D/Curve.h>
#include <vector>

using std::vector;


namespace LinearAlgebraTools {

using Math::Matrix;
using Math::Vector;
using Math::LinearSystems::solveLU;

namespace Primitives {

/**
* HermiteCurve Class
* 1D (function) Hermite Cubic curve
**/
class HermiteCurve : public Curve
{
protected:
    Vector<4> coefficients;
    bool coeff_defined;
public:
    HermiteCurve(): coeff_defined(false) {}

    HermiteCurve(const Vector<4>& coef): coefficients(coef), coeff_defined(true) {}

    virtual vector<double> coeff()
    {
        vector<double> tmp;
        tmp.insert(tmp.begin(), coefficients.data(), coefficients.data()+4);
        return tmp;
    }

protected:
    void calculateCoefficients()
    {
        if(!defined())
            return;
        double p0 = points[0], pd0 = dot_points[0][1];
        double p1 = points[1], pd1 = dot_points[1][1];
        Matrix<4,4> A{0.0, 0.0, 0.0, 1.0,
                      1.0, 1.0, 1.0, 1.0,
                      0.0, 0.0, 1.0, 0.0,
                      3.0, 2.0, 1.0, 0.0};
        if(p0!=p1)
        {
            pd0 /= (p1-p0);
            pd1 /= (p1-p0);
        }
        Vector<4> y{p0, p1, pd0, pd1};
        coefficients = solveLU(A, y);
    }

    bool defined() const
    {
        if(coeff_defined)
            return true;
        return (points.size()==2 && dot_points.size()==2);
    }

    bool canAddPoint(const double& point)
    {
        return false;
    }

    bool canAddDotPoint(const Vector<2>& point)
    {
        return true;
    }

    bool canAddDDotPoint(const Vector<2>& point)
    {
        return false;
    }
};

} }

#endif
