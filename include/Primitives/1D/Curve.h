#ifndef CURVE_H
#define CURVE_H

/**
* Includes
**/
#include <Math/Vector.h>
#include <Misc/Helper.h>
#include <vector>
#include <algorithm>

using std::vector;
using std::find;


namespace LinearAlgebraTools {

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
    vector<double> points;
    vector<Vector<2> > dot_points;
    vector<Vector<2> > ddot_points;
protected:
    virtual void calculateCoefficients() = 0;

    virtual bool defined() const = 0;

    virtual bool canAddPoint(const double& point) = 0;

    virtual bool canAddDotPoint(const Vector<2>& point) = 0;

    virtual bool canAddDDotPoint(const Vector<2>& point) = 0;
public:
    void addPoint(const double& point)
    {
        if(defined() || find(points.begin(), points.end(), point) != points.end() || !canAddPoint(point))
            return;
        points.push_back(point);
        if(defined())
            calculateCoefficients();
    }

    void addDotPoint(const Vector<2>& point)
    {
        if(defined() || find(dot_points.begin(), dot_points.end(), point) != dot_points.end() || !canAddDotPoint(point))
            return;
        points.push_back(point[0]);
        dot_points.push_back(point);
        if(defined())
            calculateCoefficients();
    }

    void addDDotPoint(const Vector<2>& point)
    {
        if(defined() || find(ddot_points.begin(), ddot_points.end(), point) != ddot_points.end() || !canAddDDotPoint(point))
            return;
        points.push_back(point[0]);
        ddot_points.push_back(point);
        if(defined())
            calculateCoefficients();
    }

    virtual vector<double> coeff() = 0;

    vector<double>& getPoints() { return points; }
    vector<Vector<2> >& getDotPoints() { return dot_points; }
    vector<Vector<2> >& getDDotPoints() { return ddot_points; }
};

} }

#endif
