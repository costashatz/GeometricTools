#ifndef CARDINAL_CURVE_H
#define CARDINAL_CURVE_H

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
* CardinalCurve Class
* 1D (function) Cardinal Cubic curve
**/
class CardinalCurve : public Curve
{
protected:
    Vector<4> coefficients;
    bool coeff_defined;
    double t_;
public:
    CardinalCurve(const double& t = 0.0): t_(t), coeff_defined(false) {}

    CardinalCurve(const Vector<4>& coef, const double& t = 0.0): t_(t), coefficients(coef), coeff_defined(true) {}

    double& t()
    {
        return t_;
    }

    void setT(const double& t)
    {
        t_ = t;
    }

    virtual double getPoint(const double &u)
    {
        double s = 0.0;
        for(int i=0;i<=3;i++)
        {
            s += coefficients[i]*std::pow(u, double(3-i));
        }
        return s;
    }

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
        double p0 = points[1], pd0 = 0.5*(1.0-t_)*(points[2]-points[0]);
        double p1 = points[2], pd1 = 0.5*(1.0-t_)*(points[3]-points[1]);
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
        return (points.size()==4);
    }

    bool canAddPoint(const double& point)
    {
        return true;
    }

    bool canAddDotPoint(const Vector<2>& point)
    {
        return false;
    }

    bool canAddDDotPoint(const Vector<2>& point)
    {
        return false;
    }
};

} }

#endif
