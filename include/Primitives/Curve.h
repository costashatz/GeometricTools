#ifndef CURVE_H
#define CURVE_H

/**
* Includes
**/
#include <Math/Vector.h>
#include <Math/Matrix.h>
#include <vector>

using std::vector;


namespace LinearAlgebraTools {

using Math::Matrix;
using Math::Vector;

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

    double mapXtoU(const double& x)
    {
        double xx = std::min(max_x, std::max(min_x, x));
        return min_u+max_u*(xx-min_x)/(max_x-min_x);
    }

    double mapUtoX(const double& u)
    {
        double uu = std::min(max_u, std::max(min_u, u));
        return min_x+max_x*(uu-min_u)/(max_u-min_u);
    }
public:
    virtual void addPoint(const Vector<2>& point)
    {
        addPointU({mapXtoU(point[0]), point[1]});
    }

    virtual void addPointU(const Vector<2>& point) = 0;

    virtual double getPoint(const double& x)
    {
        return getPointU(mapXtoU(x));
    }

    virtual double getPointU(const double& u) = 0;

    virtual void addDotPoint(const Vector<2>& point)
    {
        addDotPointU({mapXtoU(point[0]), point[1]});
    }

    virtual void addDotPointU(const Vector<2>& point) = 0;

    virtual double getDotPoint(const double& x)
    {
        return getDotPointU(mapXtoU(x));
    }

    virtual double getDotPointU(const double& x) = 0;

    virtual void addDDotPoint(const Vector<2>& point)
    {
        addDDotPointU({mapXtoU(point[0]), point[1]});
    }

    virtual void addDDotPointU(const Vector<2>& point) = 0;

    virtual double getDDotPoint(const double& x)
    {
        return getDDotPointU(mapXtoU(x));
    }

    virtual double getDDotPointU(const double& x) = 0;

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
