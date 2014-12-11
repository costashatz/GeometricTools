#ifndef PIECEWISE_CURVE_H
#define PIECEWISE_CURVE_H

/**
* Includes
**/
#include <Primitives/1D/Curve.h>

#include <Math/Vector.h>
#include <Math/Matrix.h>
#include <vector>

using std::vector;


namespace LinearAlgebraTools {

using Math::Matrix;
using Math::Vector;

namespace Primitives {

/**
* PiecewiseCurve Class
* Abstract 1D piecewise curve - base class for all 1D piecewise curves
**/
class PiecewiseCurve : public Curve
{
protected:
    vector<Curve*> curves;
public:
    PiecewiseCurve() {}

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
