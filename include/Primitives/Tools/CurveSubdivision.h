#ifndef CURVE_SUBDIVISION_H
#define CURVE_SUBDIVISION_H

/**
* Includes
**/
#include <Primitives/2D/PlaneCurve.h>
#include <Primitives/2D/Polyline.h>

namespace LinearAlgebraTools { namespace Primitives {

template<class CurveA, class CurveB>
Polyline uniformSubdivision(PlaneCurve<CurveA,CurveB>& curve, const int& N)
{
    Polyline res;
    for (int i=0;i<=N;i++)
    {
        double u = double(i)/double(N); // assumes curves goes in [0,1]
        res.addPoint(curve.getPoint(u));
    }
    return res;
}

} }

#endif
