#ifndef CURVE_H
#define CURVE_H

/**
* Includes
**/
#include <Math/Vector.h>
#include <Math/Matrix.h>
#include <Misc/Helper.h>
#include <vector>
#include <algorithm>

using std::vector;
using std::find;


namespace LinearAlgebraTools {

using Math::Matrix;
using Math::Vector;
using namespace Helper;

namespace Primitives {

/**
* Curve Class
* Generic 1D (function) curve - base class for all parametric curves
**/
class Curve
{
protected:
    vector<Vector<2> > points;
    vector<Vector<2> > dot_points;
    vector<Vector<2> > ddot_points;
    double min_u;
    double max_u;
    double min_x;
    double max_x;
protected:
    virtual void calculateCoefficients() = 0;

    virtual bool defined() const = 0;

    virtual bool canAddPoint(const Vector<2>& point) = 0;

    virtual bool canAddDotPoint(const Vector<2>& point) = 0;

    virtual bool canAddDDotPoint(const Vector<2>& point) = 0;

    double mapXtoU(const double& x)
    {
        double xx = clamp(x, min_x, max_x);
        return min_u+max_u*(xx-min_x)/(max_x-min_x);
    }

    double mapUtoX(const double& u)
    {
        double uu = clamp(u, min_u, max_u);
        return min_x+max_x*(uu-min_u)/(max_u-min_u);
    }
public:
    void addPoint(const Vector<2>& point)
    {
        addPointU({mapXtoU(point[0]), point[1]});
    }

    void addPointU(const Vector<2>& point)
    {
        if(defined() || find(points.begin(), points.end(), point) != points.end() || !canAddPoint(point))
            return;
        double u = clamp(point[0], min_u, max_u);
        points.push_back({u, point[1]});
        if(defined())
            calculateCoefficients();
    }

    double getPoint(const double& x)
    {
        return getPointU(mapXtoU(x));
    }

    virtual double getPointU(const double& u) = 0;

    void addDotPoint(const Vector<2>& point)
    {
        addDotPointU({mapXtoU(point[0]), point[1]});
    }

    void addDotPointU(const Vector<2>& point)
    {
        if(defined() || find(dot_points.begin(), dot_points.end(), point) != dot_points.end() || !canAddDotPoint(point))
            return;
        double u = clamp(point[0], min_u, max_u);
        dot_points.push_back({u, point[1]});
        if(defined())
            calculateCoefficients();
    }

    double getDotPoint(const double& x)
    {
        return getDotPointU(mapXtoU(x));
    }

    virtual double getDotPointU(const double& u) = 0;

    void addDDotPoint(const Vector<2>& point)
    {
        addDDotPointU({mapXtoU(point[0]), point[1]});
    }

    void addDDotPointU(const Vector<2>& point)
    {
        if(defined() || find(ddot_points.begin(), ddot_points.end(), point) != ddot_points.end()  || !canAddDDotPoint(point))
            return;
        double u = clamp(point[0], min_u, max_u);
        ddot_points.push_back({u, point[1]});
        if(defined())
            calculateCoefficients();
    }

    double getDDotPoint(const double& x)
    {
        return getDDotPointU(mapXtoU(x));
    }

    virtual double getDDotPointU(const double& u) = 0;

    virtual vector<double> coeff() = 0;

    double& maxU() { return max_u; }
    double& minU() { return min_u; }

    double& maxX() { return max_x; }
    double& minX() { return min_x; }

    vector<Vector<2> >& getPoints() { return points; }
    vector<Vector<2> >& getDotPoints() { return dot_points; }
    vector<Vector<2> >& getDDotPoints() { return ddot_points; }
};

} }

#endif
