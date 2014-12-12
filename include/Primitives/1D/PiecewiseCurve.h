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
        if(curves.size()==0)
            return vector<double>();
        vector<double> coef = curves[0]->coeff();
        for(int i=1;i<curves.size();i++)
        {
            auto tmp = curves[i]->coeff();
            coef.insert(coef.end(), tmp.begin(), tmp.end());
        }
        return coef;
    }

    double getPoint(const double &u)
    {
        // Assumes that all the piecewise curve is in u[0,1]
        int c = getCurveFromU(u*curves.size());
        if(c==-1)
            return 0.0;
        return curves[c]->getPoint(u*curves.size()-c);
    }
protected:
    int getCurveFromU(const double& u)
    {
        if(curves.size()==0)
            return -1;
        int c = int(u);
        if(c>=curves.size())
            c = curves.size()-1;
        return c;
    }
};

} }


#endif
