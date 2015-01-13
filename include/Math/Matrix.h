/**
Copyright (c) 2014, Konstantinos Chatzilygeroudis
All rights reserved.
Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer
    in the documentation and/or other materials provided with the distribution.
3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived
    from this software without specific prior written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING,
BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
**/

#ifndef GEOMETRIC_TOOLS_MATH_MATRIX_H
#define GEOMETRIC_TOOLS_MATH_MATRIX_H

/**
* Includes
**/
#include <iostream>
#include <cassert>
#include <Math/Vector.h>

extern "C" int dgetrf_(int *m, int *n, double *a, int * lda, int *ipiv, int *info);
extern "C" int dgetri_(int *n, double *a, int *lda, int *ipiv, double *work, int *lwork, int *info);

namespace GeometricTools { namespace Math {

/**
* General Matrix Class
* Supports NxM-dimension matrices
**/
template<unsigned int ROWS, unsigned int COLS>
class Matrix
{
protected:
    // matrix data
    double values_[COLS*ROWS+1];
    // helper variable for initialization
    unsigned int index_;

    /**
    * Initialize Matrix from double values
    * C++11 feature - variadic templates for multiple arguments
    **/
    template<typename... Args>
    void initialize(double h, Args&&... args)
    {
        if(index_>=(COLS*ROWS))
            index_=0;
        // add first argument to the matrix data
        values_[index_] = h;
        index_++;
        // Initialize with rest of arguments
        initialize(std::forward<Args>(args)...);
    }

    /**
    * Initialize Matrix from double values
    * C++11 feature - variadic templates for multiple arguments
    * Terminating Case - end of arguments
    **/
    template<typename... Args>
    void initialize()
    {
        index_ = 0;
    }

public:	

    /**
    * Default Constructor
    * initialize all data to zero
    **/
    Matrix()
    {
        memset(values_,0,COLS*ROWS*sizeof(double));
    }

    /**
    * Copy Constructor
    * @param other - Matrix to copy from
    **/
    Matrix(const Matrix& other)
    {
        memcpy(values_,other.values_,COLS*ROWS*sizeof(double));
    }

    /**
    * Constructor
    * @param h - list of double values
    **/
    template<typename... Args>
    Matrix(double h, Args&&... args)
    {
        index_ = 0;
        memset(values_,0,COLS*ROWS*sizeof(double));
        initialize(h, args...);
    }

    /**
    * Make Identity
    * Applies only to square matrices
    * If not square, nothing happens
    **/
    void identity()
    {
        if(COLS!=ROWS)
            return;
        memset(values_,0,COLS*ROWS*sizeof(double));
        for(unsigned int i=0;i<COLS;i++)
            (*this)(i,i) = 1.0;
    }

    /**
    * Make Zero
    **/
    void zero()
    {
        memset(values_,0,COLS*ROWS*sizeof(double));
    }

    /**
    * Overloading () operator
    * Access Matrix Matlab-like
    * @param i - row index to return
    * @param j - column index to return
    * @return double - value of i-th element
    **/
    double& operator()(unsigned int i, unsigned int j) //assert legal index
    {
        return values_[i*COLS+j];
    }

    /**
    * Overloading () operator - const version
    * Access Matrix Matlab-like
    * @param i - row index to return
    * @param j - column index to return
    * @return double - value of i-th element
    **/
    double operator()(unsigned int i, unsigned int j) const //assert legal index
    {
        return values_[i*COLS+j];
    }

    /**
    * Overloading *= operator
    * Multiplication of 2 Matrices
    * @param other - Matrix to multiply with
    * @return Matrix - the result
    **/
    Matrix operator*=(const Matrix<COLS,ROWS>& other)
    {
        *this = (*this)*other;
        return *this;
    }

    /**
    * Overloading += operator
    * Addition of 2 Matrices
    * @param other - Matrix to perform addition with
    * @return Matrix - the result of the addition
    **/
    Matrix operator+=(const Matrix& other)
    {
        unsigned int N = COLS*ROWS;
        cblas_daxpy(N, 1.0, &other.values_[0], 1, &values_[0], 1);
        return *this;
    }

    /**
    * Overloading -= operator
    * Subtration of 2 Matrices
    * @param other - Matrix to perform subtraction with
    * @return Matrix - the result
    **/
    Matrix operator-=(const Matrix& other)
    {
        unsigned int N = COLS*ROWS;
        cblas_daxpy(N, -1.0, &other.values_[0], 1, &values_[0], 1);
        return *this;
    }

    /**
    * Overloading *= operator
    * Multiplication with scalar (double)
    * @param other - double to multiply with
    * @return Matrix - the result
    **/
    Matrix operator*=(const double& other)
    {
        unsigned int N = COLS*ROWS;
        cblas_dscal(N, other, &values_[0], 1);
        return *this;
    }

    /**
    * Overloading /= operator
    * Division with scalar (double) - if zero ignores division (returns self)
    * @param other - double to divide with
    * @return Matrix - the result
    **/
    Matrix operator/=(const double& other)
    {
        if(std::abs(other) < std::numeric_limits<double>::epsilon())
            return (*this);
        unsigned int N = COLS*ROWS;
        cblas_dscal(N, 1.0/other, &values_[0], 1);
        return *this;
    }

    /**
    * Get Transpose of Matrix
    * @return Matrix - the transposed matrix
    **/
    Matrix<COLS,ROWS> transpose() const
    {
        Matrix<COLS,ROWS> r;
        for(unsigned int i=0;i<ROWS;i++)
        {
            for(unsigned int j=0;j<COLS;j++)
                r(j,i) = (*this)(i,j);
        }
        return r;
    }

    Matrix<COLS,ROWS> operator~() const
    {
        Matrix<COLS,ROWS> r = this->transpose();
        return r;
    }

    /**
    * Overloading + operator
    * Addition of 2 Matrices
    * @return Matrix - the result of the addition
    **/
    Matrix operator+(const Matrix<ROWS,COLS>& mat2) const
    {
        Matrix<ROWS,COLS> t = (*this);
        t += mat2;
        return t;
    }

    /**
    * Overloading - operator
    * Subtraction of 2 Matrices
    * @return Matrix - the result of the subtraction
    **/
    Matrix operator-(const Matrix<ROWS,COLS>& mat2) const
    {
        Matrix<ROWS,COLS> t = (*this);
        t -= mat2;
        return t;
    }

    /**
    * Overloading * operator
    * Multiplication with scalar (double)
    * @param other - double to multiply with
    * @return Matrix - the result
    **/
    Matrix operator*(const double& other) const
    {
        Matrix<ROWS,COLS> t = (*this);
        t *= other;
        return t;
    }

    /**
    * Overloading / operator
    * Division with scalar (double) - if zero ignores division (returns self)
    * @param other - double to divide with
    * @return Matrix - the result
    **/
    Matrix operator/(const double& other) const
    {
        Matrix<ROWS,COLS> t = (*this);
        t /= other;
        return t;
    }

    /**
    * Overloading minus (-) operator
    **/
    Matrix operator-() const
    {
        return -1*(*this);
    }

    /**
    * Overloading * operator
    * Multiplication with Vector
    **/
    Vector<COLS> operator*(const Vector<COLS>& r2) const
    {
        Vector<COLS> res;
        cblas_dgemv(CblasRowMajor, CblasNoTrans, ROWS, COLS, 1.0, values_, ROWS, r2.values_, 1, 1.0, res.values_, 1);
        return res;
    }

    /**
    * Overloading * operator
    * Multiplication with Matrix
    **/
    template<unsigned int R2>
    Matrix<ROWS,R2> operator*(const Matrix<COLS,R2>& r2)
    {
        Matrix<ROWS,R2> res;
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, ROWS, R2, COLS, 1.0, values_, ROWS, r2.values_, COLS, 0.0, res.values_, ROWS);
        return res;
    }

    /**
    * Overloading == operator
    * @param other - Matrix to compare
    * @return bool
    **/
    bool operator==(const Matrix& other) const
    {
        Matrix tmp = (*this)-other;
        return (tmp.norm()<std::numeric_limits<double>::epsilon());
    }

    /**
    * Overloading != operator
    * @param other - Matrix to compare
    * @return bool
    **/
    bool operator!=(const Matrix& other) const
    {
        Matrix tmp = (*this)-other;
        return (tmp.norm()>=std::numeric_limits<double>::epsilon());
    }

    // !-- THIS NEEDS TO BE REMOVED --! //
    /**
    * Get pointer to values array
    * @return pointer to array
    **/
    double* data()
    {
        return values_;
    }

    /**
    * Get Norm of Matrix
    * @return double - the L2 Norm of the Matrix
    **/
    double norm()
    {
        unsigned int N = COLS*ROWS;
        return cblas_dnrm2(N, &values_[0], 1);
    }

    /**
    * Get ith Row
    * @param i - row index i
    **/
    Vector<COLS> getRow(const unsigned int& i)
    {
        Vector<COLS> tmp;
        memcpy(tmp.data(), &values_[i*COLS], COLS*sizeof(double));
        return tmp;
    }

    /**
    * Set ith Row
    * @param vec - Vector
    **/
    void setRow(const unsigned int& i, const Vector<COLS>& vec)
    {
        memcpy(&values_[i*COLS], vec.data(), COLS*sizeof(double));
    }

    /**
    * Get ith Column
    * @param i - column index i
    **/
    Vector<ROWS> getCol(const unsigned int& i)
    {
        Vector<COLS> tmp;
        for(unsigned int j=0;j<ROWS;j++)
            tmp[j] = values_[j*COLS+i];
        return tmp;
    }

    /**
    * Set ith Column
    * @param vec - Vector
    **/
    void setCol(const unsigned int& i, const Vector<ROWS>& vec)
    {
        for(unsigned int j=0;j<ROWS;j++)
            values_[j*COLS+i] = vec[j];
    }

    /**
    * Swap Row i with Row j
    * @param i - row index i
    * @param j - row index j
    **/
    void swapRows(unsigned int i, unsigned int j)
    {
        Vector<COLS> t;
        for(unsigned k=0;k<COLS;k++)
            t(k) = (*this)(i,k);
        for(unsigned k=0;k<COLS;k++)
        {
            (*this)(i,k) = (*this)(j,k);
            (*this)(j,k) = t(k);
        }
    }

    /**
    * Swap Column i with Column j
    * @param i - column index i
    * @param j - column index j
    **/
    void swapCols(unsigned int i, unsigned int j)
    {
        Vector<ROWS> t;
        for(unsigned k=0;k<ROWS;k++)
            t(k) = (*this)(k,i);
        for(unsigned k=0;k<COLS;k++)
        {
            (*this)(k,i) = (*this)(k,j);
            (*this)(k,j) = t(k);
        }
    }

    /**
    * Get # rows of Matrix
    * @return unsigned int - # rows
    **/
    unsigned int rows()
    {
        return ROWS;
    }

    /**
    * Get # columns of Matrix
    * @return unsigned int - # columns
    **/
    unsigned int cols()
    {
        return COLS;
    }
};

/**
* Overloading << operator
* "print" matrix to stream
* @param ostream - stream to print the matrix to
* @param Matrix - matrix to print
**/
template <unsigned int ROWS, unsigned int COLS>
std::ostream& operator<<(std::ostream& os, const Matrix<ROWS, COLS>& obj)
{
    for(unsigned int i=0;i<ROWS-1;i++)
    {
        for(unsigned int j=0;j<COLS;j++)
            os<<obj(i,j)<<" ";
        os<<std::endl;
    }
    for(unsigned int j=0;j<COLS;j++)
        os<<obj(ROWS-1,j)<<" ";
    return os;
}

/**
* Overloading >> operator
* read matrix from stream
* @param istream - stream to read the matrix from
* @param Matrix - matrix to read to
**/
template <unsigned int ROWS, unsigned int COLS>
std::istream& operator>>(std::istream& in, Matrix<ROWS, COLS>& obj)
{
    for(unsigned int i=0;i<ROWS;i++)
    {
        for(unsigned int j=0;j<COLS;j++)
            in>>obj(i,j);
    }
    return in;
}

/**
* Overloading / operator
* Division with scalar (double)
* Perform Matrix / scalar (as opposed to Matrix / scalar inside Class)
**/
template<unsigned int N, unsigned int M>
Matrix<N,M> operator/(const double& a, const Matrix<N,M>& b)
{
    if(std::abs(a) < std::numeric_limits<double>::epsilon())
        return Matrix<N,M>();
    Matrix<N,M> tmp = Matrix<N,M>(b);
    unsigned int S = N*M;
    for(unsigned int i=0;i<S;i++) {
        if(std::abs(tmp.data()[i]) < std::numeric_limits<double>::epsilon())
            return b;
        tmp.data()[i] = a/tmp.data()[i];
    }
    return tmp;
}

/**
* Overloading * operator
* Multiplication with scalar (double)
* Perform scalar * Matrix (as opposed to Matrix * scalar inside Class)
**/
template<unsigned int N, unsigned int M>
Matrix<N,M> operator*(const double& a, Matrix<N,M>& b)
{
    return b*a;
}

/**
* Overloading * operator
* Multiplication with Vector - (vT*M)T
**/
template<unsigned int C1, unsigned int K>
Vector<C1> operator*(const Vector<K>& v1, const Matrix<K, C1>& v2)
{
    Vector<C1> res;
    cblas_dgemv(CblasRowMajor, CblasTrans, C1, K, 1.0, v1.data(), C1, v2.data(), 1, 1.0, res.data(), 1);
    return res;
}

/**
* Get Minor of the Matrix
**/
template<unsigned int ROWS, unsigned int COLS>
Matrix<ROWS-1,COLS-1> Minor(Matrix<ROWS,COLS> m, unsigned int i, unsigned int j)
{
    Matrix<ROWS-1,COLS-1> res;
    unsigned int rk=0, ck=0;
    for(unsigned int r=0;r<ROWS;r++)
    {
        if(r==i)
            continue;
        ck = 0;
        for(unsigned int c=0;c<COLS;c++)
        {
            if(c==j)
                continue;
            res(rk,ck) = m(r,c);
            ck++;
        }
        rk++;
    }
    return res;
}

/**
* Get Inverse of Matrix
* Applies only to square matrices
**/
template<unsigned int D>
Matrix<D,D> inverse(const Matrix<D,D>& mat)
{
    Matrix<D,D> res = mat;
    int d = D;
    int* ipiv = new int[D], info;
    double* work = new double[D];
    dgetrf_(&d, &d, res.data(), &d, ipiv, &info);
    dgetri_(&d, res.data(), &d, ipiv, work, &d, &info);
    delete ipiv;
    return res;
}

/**
* Get Determinant of the Matrix
**/
template<unsigned int ROWS, unsigned int COLS>
double determinant(Matrix<ROWS,COLS> m)
{
    if(COLS!=ROWS || COLS<2)
        return 0.0;
    if(COLS==2)
    {
        double x = (m(0,0)*m(1,1)-m(1,0)*m(0,1));
        return x;
    }
    double s = 0.0;
    const int R = ROWS-1;
    const int C = COLS-1;
    Matrix<R,C> temp;
    for(unsigned int i=0;i<COLS;i++)
    {
        temp = Minor(m,0,i);
        s += pow(-1,i)*m(0,i)*determinant(temp);
    }
    return s;
}

/**
* Get Determinant of the Matrix
* Base Case - 1x1 Matrix - Scalar value
**/
template<>
double determinant(Matrix<1,1> m)
{
    return 0.0;
}

/**
* Overloading ~ operator
* Vector tranpose - Returns Matrix
**/
template<unsigned int N>
Matrix<N,1> operator~(const Vector<N>& vec)
{
    Matrix<N,1> res;
    for(unsigned int i=0;i<N;i++)
        res(i,0) = vec[i];
    return res;
}

/**
* Overloading * operator
* Multiplication Vector with Vector Transpose
**/
template<unsigned int N>
Matrix<N,N> operator*(const Vector<N>& v1, const Matrix<N,1>& v2)
{
    Matrix<N,N> res;
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, N, 1, 1.0, v1.data(), 1, v2.data(), N, 0.0, res.data(), N);
    return res;
}

} }

#endif
