#ifndef SPACE_CURVE_H
#define SPACE_CURVE_H

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

template <unsigned int N>
class PolynomialCurve;

template <class CurveTypeA, class CurveTypeB>
class SpaceCurve
{
protected:
    CurveTypeA* curve_x;
    CurveTypeB* curve_y;
    double min_u;
    double max_u;
public:
    SpaceCurve(const double& minU=0.0, const double maxU=1.0): min_u(minU), max_u(maxU)
    {
        curve_x = new CurveTypeA(minU, maxU);
        curve_y = new CurveTypeB(minU, maxU);
    }

    void addPoint(const Vector<3>& point)
    {
        double u = clamp(point[0], min_u, max_u);
        curve_x->addPoint({u, point[1]});
        curve_y->addPoint({u, point[2]});
    }

    Vector<2> getPoint(const double& u)
    {
        double uu = clamp(u, min_u, max_u);
        return {curve_x->getPoint(uu), curve_y->getPoint(uu)};
    }

    void addDotPoint(const Vector<3>& point)
    {
        double u = clamp(point[0], min_u, max_u);
        curve_x->addDotPoint({u, point[1]});
        curve_y->addDotPoint({u, point[2]});
    }

    Vector<2> getDotPoint(const double& u)
    {
        double uu = clamp(u, min_u, max_u);
        return {curve_x->getDotPoint(uu), curve_y->getDotPoint(uu)};
    }

    void addDDotPoint(const Vector<3>& point)
    {
        double u = clamp(point[0], min_u, max_u);
        curve_x->addDDotPoint({u, point[1]});
        curve_y->addDDotPoint({u, point[2]});
    }

    Vector<2> getDDotPoint(const double& u)
    {
        double uu = clamp(u, min_u, max_u);
        return {curve_x->getDDotPoint(uu), curve_y->getDDotPoint(uu)};
    }

    vector<double> coeff(unsigned int i)
    {
        if(i==0)
            return curve_x->coeff();
        return curve_y->coeff();
    }

    double& maxU() { return max_u; }
    double& minU() { return min_u; }
};

/**
* Typedefs for frequently used types
**/
typedef SpaceCurve<PolynomialCurve<0>,PolynomialCurve<0> > ConstantSpaceCurve;

typedef SpaceCurve<PolynomialCurve<1>,PolynomialCurve<1> > LinearSpaceCurve;

typedef SpaceCurve<PolynomialCurve<2>,PolynomialCurve<2> > QuadraticSpaceCurve;

typedef SpaceCurve<PolynomialCurve<3>,PolynomialCurve<3> > CubicSpaceCurve;

} }

#endif
