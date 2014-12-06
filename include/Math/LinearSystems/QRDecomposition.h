#ifndef QR_DECOMPOSITION
#define QR_DECOMPOSITION

/**
* Includes
**/
#include <Math/Matrix.h>

namespace LinearAlgebraTools {

namespace Math {

namespace LinearSystems {

/**
* Perform QR Decomposition
* @param a - Matrix to decompose
* @param Q - Q Matrix --
* @param R - R Matrix  |
*					   --> Pass by reference (returns)
**/
template<unsigned int D>
void QRDecomposition(const Matrix<D,D>& a, Matrix<D,D>& Q, Matrix<D,D>& R)
{
    Q.identity();
    R.identity();
    Vector<D>* u = new Vector<D>[D];
    Vector<D>* e = new Vector<D>[D];
    for(unsigned int i=0;i<D;i++)
    {
        u[i] = a.getCol(i);
        for(unsigned int j=0;j<i;j++)
            u[i] -= projection(a.getCol(i), e[j]);
        e[i] = u[i].normalized();
        Q.setCol(i, e[i]);
    }
    for(unsigned int i=0;i<D;i++)
    {
        for(unsigned int j=i;j<D;j++)
        {
            R(i,j) = e[i]*a.getCol(j);
        }
    }
}

} } }

#endif
