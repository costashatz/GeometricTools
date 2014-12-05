#ifndef LU_DECOMPOSITION
#define LU_DECOMPOSITION
#include <Math/Matrix.h>

namespace LinearAlgebraTools {

namespace Math {

template<unsigned int ROWS, unsigned int COLS>
class Matrix;

namespace LinearSystems {

/**
* Perform LU Decomposition
* @param a - Matrix to decompose
* @param L - L Matrix --
* @param U - U Matrix  |
* @param P - P Matrix  |
*					   --> Pass by reference (returns)
**/
template<unsigned int D>
void LUDecomposition(const Matrix<D,D>& a, Matrix<D,D>& L, Matrix<D,D>& U, Matrix<D,D>& P)
{
    L.identity();
    P.identity();
    U = a;
    for(int k=0;k<D-1;k++)
    {
        double m = U(k,k);
        int newJ = k;
        for(int i=k;i<D;i++)
        {
            if(abs(U(i,k))>m)
            {
                m = abs(U(i,k));
                newJ = i;
            }
        }
        if(newJ!=k)
        {
            for(int j=k;j<D;j++)
            {
                double t = U(k,j);
                U(k,j) = U(newJ,j);
                U(newJ,j) = t;
            }
            for(int j=0;j<=k-1;j++)
            {
                double t = L(k,j);
                L(k,j) = L(newJ,j);
                L(newJ,j) = t;
            }
            P.swapRows(k,newJ);
        }

        for(int j=k+1;j<D;j++)
        {
            L(j,k) = U(j,k)/U(k,k);
            for(int i=k;i<D;i++)
            {
                U(j,i) -= L(j,k)*U(k,i);
            }
        }
    }
}

} } }

#endif
