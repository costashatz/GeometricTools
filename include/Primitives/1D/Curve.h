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
protected:
    virtual void calculateCoefficients() = 0;

    virtual bool defined() const = 0;

    virtual bool canAddPoint(const Vector<2>& point) = 0;

    virtual bool canAddDotPoint(const Vector<2>& point) = 0;

    virtual bool canAddDDotPoint(const Vector<2>& point) = 0;
public:
    void addPoint(const Vector<2>& point)
    {
        if(defined() || find(points.begin(), points.end(), point) != points.end() || !canAddPoint(point))
            return;
        points.push_back(point);
        if(defined())
            calculateCoefficients();
    }

    virtual double getPoint(const double& u) = 0;

    void addDotPoint(const Vector<2>& point)
    {
        if(defined() || find(dot_points.begin(), dot_points.end(), point) != dot_points.end() || !canAddDotPoint(point))
            return;
        dot_points.push_back(point);
        if(defined())
            calculateCoefficients();
    }

    virtual double getDotPoint(const double& u) = 0;

    void addDDotPoint(const Vector<2>& point)
    {
        if(defined() || find(ddot_points.begin(), ddot_points.end(), point) != ddot_points.end() || !canAddDDotPoint(point))
            return;
        ddot_points.push_back(point);
        if(defined())
            calculateCoefficients();
    }

    virtual double getDDotPoint(const double& u) = 0;

    virtual vector<double> coeff() = 0;

    vector<Vector<2> >& getPoints() { return points; }
    vector<Vector<2> >& getDotPoints() { return dot_points; }
    vector<Vector<2> >& getDDotPoints() { return ddot_points; }
};

} }

#endif
