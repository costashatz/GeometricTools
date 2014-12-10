#ifndef POLYNOMIAL_CURVE_H
#define POLYNOMIAL_CURVE_H

/**
* Includes
**/
#include <Math/Vector.h>
#include <Math/Matrix.h>
#include <Math/LinearSystems/SolveLU.h>
#include <Primitives/1D/Curve.h>
#include <vector>

using std::vector;


namespace LinearAlgebraTools {

using Math::Matrix;
using Math::Vector;
using Math::LinearSystems::solveLU;

namespace Primitives {

/**
* PolynomialCurve Class
* 1D (function) Polynomial curve defined by a series of control points (and slopes, curvatures) OR coefficients
**/
template<unsigned int N>
class PolynomialCurve : public Curve
{
protected:
    Vector<N+1> coefficients;
    bool coeff_defined;
public:
    PolynomialCurve(const double& minU = 0.0, const double& maxU = 1.0)
    {
        min_u = minU;
        max_u = maxU;
        coeff_defined = false;
    }

    PolynomialCurve(const Vector<N+1>& coef, const double& minU = 0.0, const double& maxU = 1.0)
    {
        min_u = minU;
        max_u = maxU;
        coefficients = coef;
        coeff_defined = true;
    }

    virtual double getPoint(const double& u)
    {
        double s = 0.0;
        for(unsigned int i=0;i<=N;i++)
            s += coefficients[i]*std::pow(u, N-i);
        return s;
    }

    virtual double getDotPoint(const double& u)
    {
        double s = 0.0;
        for(int i=0;i<=N-1;i++)
            s += (N-i)*coefficients[i]*std::pow(u, N-i-1);
        return s;
    }

    virtual double getDDotPoint(const double& u)
    {
        double s = 0.0;
        for(int i=0;i<=N-2;i++)
            s += (N-i)*(N-i-1)*coefficients[i]*std::pow(u, N-i-2);
        return s;
    }

    virtual vector<double> coeff()
    {
        vector<double> tmp;
        tmp.insert(tmp.begin(), coefficients.data(), coefficients.data()+(N+1));
        return tmp;
    }

protected:
    void calculateCoefficients()
    {
        if(!defined())
            return;

        Matrix<N+1,N+1> A;
        Vector<N+1> y;
        int ddot = ddot_points.size();
        int dot = dot_points.size();
        for(int i=0;i<ddot;i++)
        {
            y[i] = ddot_points[i][1];
            for(int j=0;j<=N-2;j++)
            {
                A(i,j) = (N-j)*(N-j-1)*std::pow(ddot_points[i][0], N-j-2);
            }
        }
        for(int i=ddot;i<dot;i++)
        {
            y[i] = dot_points[i-ddot][1];
            for(int j=0;j<=N-1;j++)
            {
                A(i,j) = (N-j)*std::pow(dot_points[i-ddot][0], N-j-1);
            }
        }
        for(unsigned int i=dot+ddot;i<=N;i++)
        {
            y[i] = points[i-(dot+ddot)][1];
            for(unsigned int j=0;j<=N;j++)
            {
                A(i,j) = std::pow(points[i-(dot+ddot)][0], N-j);
            }
        }
        coefficients = solveLU(A, y);
    }

    bool defined() const
    {
        if(coeff_defined)
            return true;
        if (points.size() + dot_points.size() + ddot_points.size() <= N)
            return false;
        for(int i=0;i<=N;i++)
        {
            int d = dot_points.size();
            int dd = ddot_points.size();
            if(i>=N)
                d = 0;
            if(i>=((int)N-1))
                dd = 0;
            if((points.size()+d+dd)==0)
                return false;
        }
        return true;
    }

    bool canAddPoint(const Vector<2> &point)
    {
        return true;
    }

    bool canAddDotPoint(const Vector<2> &point)
    {
        return (int(N)-dot_points.size()) > 0;
    }

    bool canAddDDotPoint(const Vector<2> &point)
    {
        return (int(N-1)-ddot_points.size()) > 0;
    }
};

/**
* Typedefs for frequently used types
**/
typedef PolynomialCurve<0> ConstantPolynomialCurve;

typedef PolynomialCurve<1> LinearPolynomialCurve;

typedef PolynomialCurve<2> QuadraticPolynomialCurve;

typedef PolynomialCurve<3> CubicPolynomialCurve;

} }

#endif
