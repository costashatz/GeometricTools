#ifndef PIECEWISE_CURVE_H
#define PIECEWISE_CURVE_H

/**
* Includes
**/
#include <Primitives/PolynomialCurve.h>

#include <Math/Vector.h>
#include <Math/Matrix.h>
#include <Math/LinearSystems/SolveLU.h>
#include <vector>

using std::vector;


namespace LinearAlgebraTools {

using Math::Matrix;
using Math::Vector;

namespace Primitives {

/**
* PiecewiseCurve Class
* 1D piecewise curve defined by a series of control points and between each pair a separate curve is defined
**/
class PiecewiseCurve : public Curve
{
protected:
    vector<Curve*> curves;
public:
    PiecewiseCurve(const double& minU = 0.0, const double& maxU = 1.0, const double& minX = std::numeric_limits<double>::min(), const double& maxX = std::numeric_limits<double>::max())
    {
        min_u = minU;
        max_u = maxU;
        min_x = minX;
        max_x = maxX;
    }

    virtual vector<double> coeff()
    {
        vector<double> coef;
        for(int i=0;i<curves.size();i++)
            coef.insert(coef.end(), curves[i]->coeff().begin(), curves[i]->coeff().end());
        return coef;
    }
};

} }


#endif
