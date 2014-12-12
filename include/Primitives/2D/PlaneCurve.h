#ifndef PLANE_CURVE_H
#define PLANE_CURVE_H

/**
* Includes
**/
#include <Math/Vector.h>
#include <Math/Matrix.h>
#include <Math/LinearSystems/SolveLU.h>
#include <Misc/Helper.h>
#include <vector>

using std::vector;

namespace LinearAlgebraTools {

using Math::Matrix;
using Math::Vector;
using namespace Helper;

namespace Primitives {

template <class CurveTypeA, class CurveTypeB>
class PlaneCurve
{
protected:
    CurveTypeA* curve_x;
    CurveTypeB* curve_y;
public:
    template<typename...U>
    PlaneCurve(U... args)
    {
        curve_x = new CurveTypeA(std::forward<U>(args)...);
        curve_y = new CurveTypeB(std::forward<U>(args)...);
    }

    void addPoint(const Vector<2>& point)
    {
        curve_x->addPoint(point[0]);
        curve_y->addPoint(point[1]);
    }

    void addPointInPlace(const Vector<2>& point, const unsigned int& i)
    {
        curve_x->addPointInPlace(point[0], i);
        curve_y->addPointInPlace(point[1], i);
    }

    void addDotPoint(const Vector<4>& point)
    {
        curve_x->addDotPoint({point[0], point[2]});
        curve_y->addDotPoint({point[1], point[3]});
    }

    void addDotPointInPlace(const Vector<4>& point, const unsigned int& i)
    {
        curve_x->addDotPointInPlace({point[0], point[2]}, i);
        curve_y->addDotPointInPlace({point[1], point[3]}, i);
    }

    void addDDotPoint(const Vector<4>& point)
    {
        curve_x->addDDotPoint({point[0], point[2]});
        curve_y->addDDotPoint({point[1], point[3]});
    }

    void addDDotPointInPlace(const Vector<4>& point, const unsigned int& i)
    {
        curve_x->addDDotPointInPlace({point[0], point[2]}, i);
        curve_y->addDDotPointInPlace({point[1], point[3]}, i);
    }

    Vector<2> getPoint(const double& u)
    {
        return {curve_x->getPoint(u), curve_y->getPoint(u)};
    }

    vector<double> coeff(unsigned int i)
    {
        if(i==0)
            return curve_x->coeff();
        return curve_y->coeff();
    }
};

} }

#endif
