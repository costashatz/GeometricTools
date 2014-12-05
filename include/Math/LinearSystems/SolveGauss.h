#ifndef SOLVE_GAUSS_H
#define SOLVE_GAUSS_H
#include <Math/Matrix.h>

namespace LinearAlgebraTools {

namespace Math {

namespace LinearSystems {

/**
* Solve Linear System using Gauss Elimination
* @param A - parameter Matrix
* @param B - constant Vector
* @return Vector<D> - solution
**/
template<unsigned int D>
Vector<D> solveGauss(const Matrix<D,D>& A, const Vector<D>& B)
{
    Matrix<D,D> a = A;
    Vector<D> b = B;
    Vector<D> res;
    for(unsigned int i=0;i<D-1;i++)
    {
        int m = 0;
        double temp1 = 0.0;
        for(unsigned int j=i;j<D;j++)
        {
            double temp2 = 0.0;
            for(unsigned int k=i+1;k<D;k++)
                if(abs(a(j,k))>temp2)
                    temp2 = abs(a(j,k));
            if(abs(a(j,i))>(temp1*temp2))
            {
                temp1 = abs(a(j,i))/temp2;
                m = j;
            }
        }
        if(m!=i)
        {
            a.swapRows(i,m);
            double temp = b(i);
            b(i) = b(m);
            b(m) = temp;
        }

        if(abs(a(i,i))<std::numeric_limits<double>::epsilon())
            return res;

        for(unsigned int k=i+1;k<D;k++)
        {
            double r = a(k,i)/a(i,i);
            for(unsigned int c=0;c<D;c++)
            {
                a(k,c) -= r*a(i,c);
            }
            b(k) -= r*b(i);
        }
    }
    for(int k=D-1;k>=0;k--)
    {
        double S = 0.0;
        for(unsigned int j=k+1;j<D;j++)
        {
            S += a(k,j)*res(j);
        }
        res(k) = 1.0/a(k,k)*(b(k)-S);
    }
    return res;
}

} } }

#endif
