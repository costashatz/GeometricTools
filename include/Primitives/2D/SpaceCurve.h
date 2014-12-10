#ifndef SPACE_CURVE_H
#define SPACE_CURVE_H

/**
* Includes
**/
#include <Primitives/Curve.h>

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
protected:
    double mapXtoU(const double& x)
    {
        double xx = clamp(x, curve_x->minX(), curve_x->maxX());
        return min_u+max_u*(xx-curve_x->minX())/(curve_x->maxX()-curve_x->minX());
    }

    double mapUtoX(const double& u)
    {
        double uu = clamp(u, min_u, max_u);
        return curve_x->minX()+curve_x->maxX()*(uu-min_u)/(max_u-min_u);
    }

    double mapYtoU(const double& y)
    {
        double yy = clamp(y, curve_y->minX(), curve_y->maxX());
        return min_u+max_u*(yy-curve_y->minX())/(curve_y->maxX()-curve_y->minX());
    }

    double mapUtoY(const double& u)
    {
        double uu = clamp(u, min_u, max_u);
        return curve_y->minX()+curve_y->maxX()*(uu-min_u)/(max_u-min_u);
    }
public:
    SpaceCurve(const double& minU=0.0, const double maxU=1.0, const double& minX = std::numeric_limits<double>::min(), const double& maxX = std::numeric_limits<double>::max(), const double& minY = std::numeric_limits<double>::min(), const double& maxY = std::numeric_limits<double>::max()): min_u(minU), max_u(maxU)
    {
        curve_x = new CurveTypeA(minU, maxU, minX, maxX);
        curve_y = new CurveTypeB(minU, maxU, minY, maxY);
    }

    void addPoint(const Vector<4>& point)
    {
        addPointU({mapXtoU(point[0]), mapYtoU(point[1]), point[2], point[3]});
    }

    void addPointU(const Vector<4>& point)
    {
        double ux = clamp(point[0], min_u, max_u);
        double uy = clamp(point[1], min_u, max_u);
        curve_x->addPointU({ux, point[2]});
        curve_y->addPointU({uy, point[3]});
    }

    Vector<2> getPoint(const Vector<2>& x)
    {
        return getPointU(mapXtoU(x[0]), mapYtoU(x[1]));
    }

    Vector<2> getPointU(const Vector<2>& u)
    {
        double ux = clamp(u[0], min_u, max_u);
        double uy = clamp(u[1], min_u, max_u);
        return {curve_x->getPointU(ux), curve_y->getPointU(uy)};
    }

    void addDotPoint(const Vector<4>& point)
    {
        addDotPointU({mapXtoU(point[0]), mapYtoU(point[1]), point[2], point[3]});
    }

    void addDotPointU(const Vector<3>& point)
    {
        double ux = clamp(point[0], min_u, max_u);
        double uy = clamp(point[1], min_u, max_u);
        curve_x->addDotPointU({ux, point[2]});
        curve_y->addDotPointU({uy, point[3]});
    }

    Vector<2> getDotPoint(const Vector<2>& p)
    {
        return getDotPointU({mapXtoU(p[0]), mapYtoU(p[1])});
    }

    Vector<2> getDotPointU(const Vector<2>& u)
    {
        double ux = clamp(u[0], min_u, max_u);
        double uy = clamp(u[1], min_u, max_u);
        return {curve_x->getDotPointU(ux), curve_y->getDotPointU(uy)};
    }

    void addDDotPoint(const Vector<4>& point)
    {
        addDDotPointU({mapXtoU(point[0]), mapYtoU(point[1]), point[2], point[3]});
    }

    void addDDotPointU(const Vector<4>& point)
    {
        double ux = clamp(point[0], min_u, max_u);
        double uy = clamp(point[1], min_u, max_u);
        curve_x->addDDotPointU({ux, point[2]});
        curve_y->addDDotPointU({uy, point[3]});
    }

    Vector<2> getDDotPoint(const Vector<2>& p)
    {
        return getDDotPointU(mapXtoU(p[0]), mapYtoU(p[1]));
    }

    Vector<2> getDDotPointU(const Vector<2>& u)
    {
        double ux = clamp(u[0], min_u, max_u);
        double uy = clamp(u[1], min_u, max_u);
        return {curve_x->getDDotPointU(ux), curve_y->getDDotPointU(uy)};
    }

    vector<double> coeff(unsigned int i)
    {
        if(i==0)
            return curve_x->coeff();
        return curve_y->coeff();
    }

    double& maxU() { return max_u; }
    double& minU() { return min_u; }
    double& maxX() { return curve_x->maxX(); }
    double& minX() { return curve_x->minX(); }
    double& maxY() { return curve_y->maxX(); }
    double& minY() { return curve_y->minX(); }
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
