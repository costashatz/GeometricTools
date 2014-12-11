#ifndef NATURAL_PIECEWISE_CURVE_H
#define NATURAL_PIECEWISE_CURVE_H

/**
* Includes
**/
#include <Primitives/1D/PiecewiseCurve.h>
#include <Primitives/1D/PolynomialCurve.h>

#include <Math/Vector.h>
#include <Math/Matrix.h>
#include <vector>

using std::vector;


namespace LinearAlgebraTools {

using Math::Matrix;
using Math::Vector;

namespace Primitives {

/**
* NaturalPiecewiseCurve Class
* 1D natural piecewise curve - C0,C1,C2 constraints
**/
class NaturalCubicPiecewiseCurve : public PiecewiseCurve
{
protected:
    double dot_first;
    double dot_last;
public:
    NaturalCubicPiecewiseCurve(const double& dfirst, const double& dlast): dot_first(dfirst), dot_last(dlast), PiecewiseCurve() {}
protected:
    void calculateCoefficients()
    {
        curves.clear();
    }

    bool canAddPoint(const Vector<2> &point)
    {
        return true;
    }

    bool canAddDotPoint(const Vector<2> &point)
    {
        return false;
    }

    bool canAddDDotPoint(const Vector<2> &point)
    {
        return false;
    }

    bool defined() const
    {
        return points.size()>=2;
    }

public:
    double getPoint(const double &u)
    {
    }

    double getDotPoint(const double &u)
    {
    }

    double getDDotPoint(const double &u)
    {
    }
};

} }


#endif
