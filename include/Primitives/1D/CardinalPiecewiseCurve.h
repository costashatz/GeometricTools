#ifndef PIECEWISE_CARDINAL_CURVE_H
#define PIECEWISE_CARDINAL_CURVE_H

/**
* Includes
**/
#include <Math/Vector.h>
#include <Math/Matrix.h>
#include <Primitives/1D/CardinalCurve.h>
#include <Primitives/1D/PiecewiseCurve.h>
#include <vector>

using std::vector;


namespace LinearAlgebraTools {

using Math::Matrix;
using Math::Vector;

namespace Primitives {

/**
* CardinalPiecewiseCurve Class
* 1D (function) Piecewise Cardinal Cubic curve
**/
class CardinalPiecewiseCurve : public PiecewiseCurve
{
protected:
    double t_;
public:
    CardinalPiecewiseCurve(const double& t = 0.0): t_(t) {std::cout<<t<<std::endl;}

    double& t()
    {
        return t_;
    }

    void setT(const double& t)
    {
        t_ = t;
    }

    void addPoint(const double& point)
    {
        points.push_back(point);
        if(defined())
            calculateCoefficients();
    }

protected:
    void calculateCoefficients()
    {
        if(!defined())
            return;
        curves.push_back(new CardinalCurve(t_));
        for(int i=points.size()-4;i<points.size();i++)
        {
            ((CardinalCurve*)curves[curves.size()-1])->addPoint(points[i]);
        }
    }

    bool defined() const
    {
        return (points.size()>=4);
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
