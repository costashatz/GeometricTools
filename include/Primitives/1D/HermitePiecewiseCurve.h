#ifndef PIECEWISE_HERMITE_CURVE_H
#define PIECEWISE_HERMITE_CURVE_H

/**
* Includes
**/
#include <Math/Vector.h>
#include <Math/Matrix.h>
#include <Primitives/1D/HermiteCurve.h>
#include <Primitives/1D/PiecewiseCurve.h>
#include <vector>

using std::vector;


namespace LinearAlgebraTools {

using Math::Matrix;
using Math::Vector;

namespace Primitives {

/**
* HermiteCurve Class
* 1D (function) Piecewise Hermite Cubic curve
**/
class HermitePiecewiseCurve : public PiecewiseCurve
{
public:
    HermitePiecewiseCurve() {}

    void addDotPoint(const Vector<2> &point)
    {
        if(std::find(dot_points.begin(), dot_points.end(), point) != dot_points.end())
            return;
        points.push_back(point[0]);
        dot_points.push_back(point);
        if(defined())
            calculateCoefficients();
    }

protected:
    void calculateCoefficients()
    {
        if(!defined())
            return;
        for(int i=dot_points.size()-2;i<dot_points.size()-1;i++)
        {
            curves.push_back(new HermiteCurve());
            ((HermiteCurve*)curves[i])->addDotPoint(dot_points[i]);
            ((HermiteCurve*)curves[i])->addDotPoint(dot_points[i+1]);
        }
    }

    bool defined() const
    {
        return (points.size()>=2 && dot_points.size()>=2);
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
