#ifndef POLYNOMIAL_CURVE_H
#define POLYNOMIAL_CURVE_H

/**
* Includes
**/
#include <Math/Vector.h>
#include <Math/Matrix.h>
#include <Math/LinearSystems/SolveLU.h>
#include <vector>
#include <algorithm>

using std::vector;
using std::find;


namespace LinearAlgebraTools {

using Math::Matrix;
using Math::Vector;
using Math::LinearSystems::solveLU;

namespace Primitives {

/**
* PolynomialCurve Class
* 1D (function) Polynomial curve defined by a series of control points (and slopes, curvatures) OR coefficients
**/
template<unsigned int N, unsigned int D>
class PolynomialCurve
{
protected:
    Vector<N> coefficients;
    vector<Vector<D> > points;
    vector<Vector<D> > dot_points;
    vector<Vector<D> > ddot_points;
    double min_u;
    double max_u;
public:
    PolynomialCurve(const double& minU = 0.0, const double& maxU = 1.0): min_u(minU), max_u(maxU) {}

    PolynomialCurve(const Vector<N>& coef, const double& minU = 0.0, const double& maxU = 1.0): min_u(minU), max_u(maxU), coefficients(coef) {}

    virtual void addPoint(const Vector<D>& point)
    {
        if(defined() || find(points.begin(), points.end(), point) != points.end())
            return;
        points.push_back(point);
    }

    double getPoint(const double& u)
    {
        double uu = std::min(max_u, std::max(min_u, u));
        double s = 0.0;
        for(unsigned int i=0;i<N;i++)
            s += coefficients[i]*std::pow(uu, N-i-1);
        return s;
    }

    virtual void addDotPoint(const Vector<D>& point)
    {
        if(defined() || ((int)N-dot_points.size()) <= 1 || find(dot_points.begin(), dot_points.end(), point) != dot_points.end())
            return;
        dot_points.push_back(point);
    }

    double getDotPoint(const double& u)
    {
        double uu = std::min(max_u, std::max(min_u, u));
        double s = 0.0;
        for(int i=0;i<N-1;i++)
            s += (N-i-1)*coefficients[i]*std::pow(u, N-i-2);
        return s;
    }

    virtual void addDDotPoint(const Vector<D>& point)
    {
        if(defined() || ((int)N-ddot_points.size()) <= 1 || find(ddot_points.begin(), ddot_points.end(), point) != ddot_points.end())
            return;
        ddot_points.push_back(point);
    }

    double getDDotPoint(const double& u)
    {
        double uu = std::min(max_u, std::max(min_u, u));
        double s = 0.0;
        for(int i=0;i<N-2;i++)
            s += (N-i-1)*(N-i-2)*coefficients[i]*std::pow(u, N-i-3);
        return s;
    }
protected:
    void calculateCoefficients()
    {
        // NEEDS CHECKING FOR THE MATH
        if(!defined())
            return;

        Matrix<N,N> A;
        Vector<N> y;
        int ddot = ddot_points.size();
        int dot = dot_points.size();
        for(int i=0;i<ddot;i++)
        {
            y[i] = ddot_points[i][1];
            for(int j=0;j<N-2;j++)
            {
                A(i,j) = (N-j-1)*(N-j-2)*std::pow(ddot_points[i][0], N-j-3);
            }
        }
        for(int i=ddot;i<dot;i++)
        {
            y[i] = dot_points[i-ddot][1];
            for(int j=0;j<N-1;j++)
            {
                A(i,j) = (N-j-1)*std::pow(dot_points[i-ddot][0], N-j-2);
            }
        }
        for(unsigned int i=dot+ddot;i<N;i++)
        {
            y[i] = points[i-(dot+ddot)][1];
            for(unsigned int j=0;j<N;j++)
            {
                A(i,j) = std::pow(points[i-(dot+ddot)][0], N-j-1);
            }
        }
        coefficients = solveLU(A, y);
    }

    bool defined() const
    {
        if (points.size() + dot_points.size() + ddot_points.size() < N)
            return false;
        for(int i=0;i<N;i++)
        {
            int d = dot_points.size();
            int dd = ddot_points.size();
            if(i>=(N-1))
                d = 0;
            if(i>=((int)N-2))
                dd = 0;
            if((points.size()+d+dd)==0)
                return false;
        }
        return true;
    }
};

/**
* Typedefs for frequently used types
**/
typedef PolynomialCurve<1,2> ConstantPolynomialCurve2;
typedef PolynomialCurve<1,3> ConstantPolynomialCurve3;

typedef PolynomialCurve<2,2> LinearPolynomialCurve2;
typedef PolynomialCurve<2,3> LinearPolynomialCurve3;

typedef PolynomialCurve<3,2> QuadraticPolynomialCurve2;
typedef PolynomialCurve<3,3> QuadraticPolynomialCurve3;

typedef PolynomialCurve<4,2> CubicPolynomialCurve2;
typedef PolynomialCurve<4,3> CubicPolynomialCurve3;

} }

#endif
