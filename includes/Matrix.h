#include <iostream>
#include <cassert>
#include "Vector.h"
using namespace std;

#define MAX(a,b) a>b?a:b

void strassen(double **a, double **b, double **c, int tam);
void sum(double **a, double **b, double **result, int tam);
void subtract(double **a, double **b, double **result, int tam);
double **allocate_real_matrix(int tam, int random);
double **free_real_matrix(double **v, int tam);


template<unsigned int ROWS, unsigned int COLS>
class Matrix
{
protected:
	double values[COLS*ROWS+1];
	unsigned int index;

	template<typename... Args>
	void Initialize(double h, Args&&... args)
	{
		if(index>=(COLS*ROWS))
			index=0;
		unsigned int i = index/COLS;
		unsigned int j = index%COLS;
		values[i*COLS+j] = h;
		index++;
		Initialize(forward<Args>(args)...);
	}

	template<typename... Args>
	void Initialize()
	{
		index = 0;
	}

public:	

	Matrix()
	{
		memset(values,0,COLS*ROWS*sizeof(double));
	}

	Matrix(const Matrix& other)
	{
		memcpy(values,other.values,COLS*ROWS*sizeof(double));
	}

	template<typename... Args>
	Matrix(double h, Args&&... args)
	{
		index = 0;
		memset(values,0,COLS*ROWS*sizeof(double));
		Initialize(h, args...);
	}

	void Identity() //assert square
	{
		memset(values,0,COLS*ROWS*sizeof(double));
		for(unsigned int i=0;i<COLS;i++)
			(*this)(i,i) = 1.0;
	}

	double& operator()(unsigned int i, unsigned int j) //assert legal index
	{
		return values[i*COLS+j];
	}



	Matrix operator+(const Matrix& other)
	{
		unsigned int N = COLS*ROWS;
		Matrix t = Matrix(*this);
		for(unsigned int i=0;i<N;i++)
			t.values[i] += other.values[i];
		return t;
	}

	Matrix operator+=(const Matrix& other)
	{
		*this = *this+other;
		return *this;
	}

	Matrix operator-(const Matrix& other)
	{
		unsigned int N = COLS*ROWS;
		Matrix t = Matrix(*this);
		for(unsigned int i=0;i<N;i++)
			t.values[i] -= other.values[i];
		return t;
	}

	Matrix operator-=(const Matrix& other)
	{
		*this = *this-other;
		return *this;
	}

	Matrix operator*(const double& other)
	{
		unsigned int N = COLS*ROWS;
		Matrix t = Matrix(*this);
		for(unsigned int i=0;i<N;i++)
			t.values[i] *= other;
		return t;
	}

	Matrix operator*=(const double& other)
	{
		*this = *this*other;
		return *this;
	}

	Matrix operator/(const double& other)
	{
		unsigned int N = COLS*ROWS;
		Matrix t = Matrix(*this);
		for(unsigned int i=0;i<N;i++)
			t.values[i] /= other;
		return t;
	}

	Matrix operator/=(const double& other)
	{
		*this = *this/other;
		return *this;
	}

	Matrix<COLS,ROWS> Transpose()
	{
		Matrix<COLS,ROWS> r;
		for(unsigned int i=0;i<ROWS;i++)
		{
			for(unsigned int j=0;j<COLS;j++)
				r(j,i) = (*this)(i,j);
		}
		return r;
	}

	Matrix<COLS,ROWS> operator~()
	{
		Matrix<COLS,ROWS> r = this->Transpose();
		return r;
	}

	void Print()
	{
		for(unsigned int i=0;i<ROWS;i++)
		{
			for(unsigned int j=0;j<COLS;j++)
				cout<<(*this)(i,j)<<" ";
			cout<<endl;
		}
	}

	double Norm()
	{
		double s = 0.;
		for(unsigned int i=0;i<ROWS;i++)
		{
			for(unsigned int j=0;j<COLS;j++)
				s += (*this)(i,j)*(*this)(i,j);
		}
		return sqrt(s);
	}

	void SwapRows(unsigned int i, unsigned int j)
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

	void SwapCols(unsigned int i, unsigned int j)
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
};

/*template<>
class Matrix<1,1>
{
protected:
	double value;

public:
	double& operator()(unsigned int i, unsigned int j) //assert legal index
	{
		return value;
	}
};*/

template<unsigned int N, unsigned int M>
Matrix<N,M> operator*(const double& a, Matrix<N,M>& b)
{
	return b*a;
}

template<unsigned int N, unsigned int M>
Matrix<N,M> operator-(Matrix<N,M>& b)
{
	return -1*b;
}

template<unsigned int C1, unsigned int K>
Vector<C1> operator*(const Matrix<C1,K>& r1, const Vector<K>& r2)
{
	Matrix<K,1> b;
	Vector<K> v = r2;
	for(unsigned int i=0;i<K;i++)
		b(i,0) = v[i];
	Matrix<C1,1> c;
	c = r1*b;
	Vector<C1> res;
	for(unsigned int i=0;i<K;i++)
		res[i] = c(i,0);
	return res;
}


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


template<unsigned int ROWS, unsigned int COLS>
double Determinant(Matrix<ROWS,COLS> m)
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
		s += pow(-1,i)*m(0,i)*Determinant(temp);
	}
	return s;
}


template<>
double Determinant(Matrix<1,1> m)
{
	return 0.0;
}


template<unsigned int C1, unsigned int K, unsigned int R2>
Matrix<C1,R2> operator*(const Matrix<C1,K>& r1, const Matrix<K,R2>& r2)
{
	/*Matrix<C1,R2> r;
	auto a = r1;
	auto b = r2;
	for(unsigned int i=0;i<C1;i++)
	{
		for(unsigned int j=0;j<R2;j++)
		{
			for(unsigned int k=0;k<K;k++)
			{
				r(i,j) += a(i,k)*b(k,j);
			}
		}
	}

	return r;*/
	Matrix<C1,K> r = r1;
	Matrix<K,R2> rb = r2;
	double** a, **b, **c;
	double D1=(ceil(log(K)/log(2.))), D2=(ceil(log(C1)/log(2.))),D3=(ceil(log(R2)/log(2.)));
	unsigned int DIM = MAX(D1,MAX(D2,D3));
	DIM = (unsigned int)pow(2,DIM);
	a = allocate_real_matrix(DIM,0);
	b = allocate_real_matrix(DIM,0);
	c = allocate_real_matrix(DIM,0);
	for(unsigned int i=0;i<C1;i++)
	{
		for(unsigned int j=0;j<K;j++)
		{
			a[i][j] = r(i,j);
		}
	}
	for(unsigned int i=0;i<K;i++)
	{
		for(unsigned int j=0;j<R2;j++)
		{
			b[i][j] = rb(i,j);
		}
	}
	strassen(a,b,c,DIM);
	a = free_real_matrix(a,DIM);
	b = free_real_matrix(b,DIM);
	Matrix<C1,R2> res;
	for(unsigned int i=0;i<C1;i++)
	{
		for(unsigned int j=0;j<R2;j++)
		{
			res(i,j) = c[i][j];
		}
	}
	c = free_real_matrix(c,DIM);
	return res;
}



/*template<int C1, int K, int R2>
Matrix<C1,R2> Multiplication(const Matrix<C1,K>& r1, const Matrix<K,R2>& r2)
{
	if(C1==K && K==R2)
		return strassenMultiplication(r1,r2);
	const int D1=(const int)(ceil(log(K)/log(2.))), D2=(const int)(ceil(log(C1)/log(2.))),D3=(const int)(ceil(log(R2)/log(2.)));
	const int DIM = MAX(D1,MAX(D2,D3));
	//const int DIM = MAX(C1,MAX(K,R2));
	Matrix<C1,K> r = r1;
	Matrix<K,R2> rb = r2;
	Matrix<DIM,DIM> a, b, c;
	for(int i=0;i<C1;i++)
	{
		for(int j=0;j<K;j++)
		{
			a(i,j) = r(i,j);
		}
	}
	for(int i=0;i<K;i++)
	{
		for(int j=0;j<R2;j++)
		{
			b(i,j) = rb(i,j);
		}
	}
	c = strassenMultiplication(a,b);
	Matrix<C1,R2> res;
	for(int i=0;i<C1;i++)
	{
		for(int j=0;j<R2;j++)
		{
			res(i,j) = c(i,j);
		}
	}
	return res;
}*/


/*template<int D>
Matrix<D,D> strassenMultiplication(const Matrix<D,D>& r1, const Matrix<D,D>& r2)
{
	if(D==0)
	{
		Matrix<D,D> res;
		return res;
	}
	Matrix<D,D> c;
	Matrix<D,D> a = r1;
	Matrix<D,D> b = r2;
	const int newDim = D/2;
	if(D==1)
	{
		c(0,0) = a(0,0)*b(0,0);
		return c;
	}

	Matrix<newDim,newDim> a11, a12, a21, a22, b11, b12, b21, b22, c11, c12, c21, c22, p1, p2, p3, p4, p5, p6, p7, aResult, bResult;

	for (int i = 0; i < newDim; i++)
	{
		for (int j = 0; j < newDim; j++)
		{
			a11(i,j) = a(i,j);
			a12(i,j) = a(i,j + newDim);
			a21(i,j) = a(i + newDim,j);
			a22(i,j) = a(i + newDim,j + newDim);
 
			b11(i,j) = b(i,j);
			b12(i,j) = b(i,j + newDim);
			b21(i,j) = b(i + newDim,j);
			b22(i,j) = b(i + newDim,j + newDim);
		}
	}

	aResult = a11+a22;
	bResult = b11+b22;
	p1 = strassenMultiplication(aResult,bResult);


	aResult = a21+a22;
	p2 = strassenMultiplication(aResult,b11);


	bResult = b12-b22;
    p3 = strassenMultiplication(a11,bResult);

	bResult = b21-b11; 
    p4 = strassenMultiplication(a22,bResult);

	aResult = a11+a12;
    p5 = strassenMultiplication(aResult,b22);

	aResult = a21-a11;
    bResult = b11+b12;
    p6 = strassenMultiplication(aResult,bResult);

	aResult = a12-a22;
	bResult = b21+b22;
	p7 = strassenMultiplication(aResult,bResult);

	c12 = p3+p5;
	c21 = p2+p4;
	c11 = p1 + p4 - p5 + p7;
	c22 = p1 + p3 - p2 + p6;

	for (int i = 0; i < newDim ; i++)
	{
		for (int j = 0 ; j < newDim ; j++)
		{
			c(i,j) = c11(i,j);
			c(i,j + newDim) = c12(i,j);
			c(i + newDim,j) = c21(i,j);
			c(i + newDim,j + newDim) = c22(i,j);
		}
	}

	return c;
}*/



/*template<int C1, int K, int R2>
Matrix<C1,R2> Multiplication(const Matrix<C1,K>& r1, const Matrix<K,R2>& r2)
{
	Matrix<C1,K> r = r1;
	Matrix<K,R2> rb = r2;
	double** a, **b, **c;
	int D1=(ceil(log(K)/log(2.))), D2=(ceil(log(C1)/log(2.))),D3=(ceil(log(R2)/log(2.)));
	int DIM = MAX(D1,MAX(D2,D3));
	DIM = (int)pow(DIM,2);
	a = allocate_real_matrix(DIM,0);
	b = allocate_real_matrix(DIM,0);
	c = allocate_real_matrix(DIM,0);
	for(int i=0;i<C1;i++)
	{
		for(int j=0;j<K;j++)
		{
			a[i][j] = r(i,j);
		}
	}
	for(int i=0;i<K;i++)
	{
		for(int j=0;j<R2;j++)
		{
			b[i][j] = rb(i,j);
		}
	}
	strassen(a,b,c,DIM);
	a = free_real_matrix(a,DIM);
	b = free_real_matrix(b,DIM);
	Matrix<C1,R2> res;
	for(int i=0;i<C1;i++)
	{
		for(int j=0;j<R2;j++)
		{
			res(i,j) = c[i][j];
		}
	}
	c = free_real_matrix(c,DIM);
	return res;
}*/



/*void strassen(double **a, double **b, double **c, int tam);
void sum(double **a, double **b, double **result, int tam);
void subtract(double **a, double **b, double **result, int tam);
double **allocate_real_matrix(int tam, int random);
double **free_real_matrix(double **v, int tam);*/
 
void strassen(double **a, double **b, double **c, int tam) {
 
    // trivial case: when the matrix is 32 X 32 or less:
    if (tam <= 32) {
		for(int i=0;i<tam;i++)
		{
			for(int j=0;j<tam;j++)
			{
				c[i][j] = 0;
				for(int k=0;k<tam;k++)
				{
					c[i][j] += a[i][k]*b[k][j];
				}
			}
		}
        return;
    }
 
    // other cases are treated here:
        int newTam = tam/2;
        double **a11, **a12, **a21, **a22;
        double **b11, **b12, **b21, **b22;
        double **c11, **c12, **c21, **c22;
        double **p1, **p2, **p3, **p4, **p5, **p6, **p7;
 
        // memory allocation:
        a11 = allocate_real_matrix(newTam, -1);
        a12 = allocate_real_matrix(newTam, -1);
        a21 = allocate_real_matrix(newTam, -1);
        a22 = allocate_real_matrix(newTam, -1);
 
        b11 = allocate_real_matrix(newTam, -1);
        b12 = allocate_real_matrix(newTam, -1);
        b21 = allocate_real_matrix(newTam, -1);
        b22 = allocate_real_matrix(newTam, -1);
 
        c11 = allocate_real_matrix(newTam, -1);
        c12 = allocate_real_matrix(newTam, -1);
        c21 = allocate_real_matrix(newTam, -1);
        c22 = allocate_real_matrix(newTam, -1);
 
        p1 = allocate_real_matrix(newTam, -1);
        p2 = allocate_real_matrix(newTam, -1);
        p3 = allocate_real_matrix(newTam, -1);
        p4 = allocate_real_matrix(newTam, -1);
        p5 = allocate_real_matrix(newTam, -1);
        p6 = allocate_real_matrix(newTam, -1);
        p7 = allocate_real_matrix(newTam, -1);
 
        double **aResult = allocate_real_matrix(newTam, -1);
        double **bResult = allocate_real_matrix(newTam, -1);
 
        int i, j;
 
        //dividing the matrices in 4 sub-matrices:
        for (i = 0; i < newTam; i++) {
            for (j = 0; j < newTam; j++) {
                a11[i][j] = a[i][j];
                a12[i][j] = a[i][j + newTam];
                a21[i][j] = a[i + newTam][j];
                a22[i][j] = a[i + newTam][j + newTam];
 
                b11[i][j] = b[i][j];
                b12[i][j] = b[i][j + newTam];
                b21[i][j] = b[i + newTam][j];
                b22[i][j] = b[i + newTam][j + newTam];
            }
        }
 
        // Calculating p1 to p7:
 
        sum(a11, a22, aResult, newTam); // a11 + a22
        sum(b11, b22, bResult, newTam); // b11 + b22
        strassen(aResult, bResult, p1, newTam); // p1 = (a11+a22) * (b11+b22)
 
        sum(a21, a22, aResult, newTam); // a21 + a22
        strassen(aResult, b11, p2, newTam); // p2 = (a21+a22) * (b11)
 
        subtract(b12, b22, bResult, newTam); // b12 - b22
        strassen(a11, bResult, p3, newTam); // p3 = (a11) * (b12 - b22)
 
        subtract(b21, b11, bResult, newTam); // b21 - b11
        strassen(a22, bResult, p4, newTam); // p4 = (a22) * (b21 - b11)
 
        sum(a11, a12, aResult, newTam); // a11 + a12
        strassen(aResult, b22, p5, newTam); // p5 = (a11+a12) * (b22)   
 
        subtract(a21, a11, aResult, newTam); // a21 - a11
        sum(b11, b12, bResult, newTam); // b11 + b12
        strassen(aResult, bResult, p6, newTam); // p6 = (a21-a11) * (b11+b12)
 
        subtract(a12, a22, aResult, newTam); // a12 - a22
        sum(b21, b22, bResult, newTam); // b21 + b22
        strassen(aResult, bResult, p7, newTam); // p7 = (a12-a22) * (b21+b22)
 
        // calculating c21, c21, c11 e c22:
 
        sum(p3, p5, c12, newTam); // c12 = p3 + p5
        sum(p2, p4, c21, newTam); // c21 = p2 + p4
 
        sum(p1, p4, aResult, newTam); // p1 + p4
        sum(aResult, p7, bResult, newTam); // p1 + p4 + p7
        subtract(bResult, p5, c11, newTam); // c11 = p1 + p4 - p5 + p7
 
        sum(p1, p3, aResult, newTam); // p1 + p3
        sum(aResult, p6, bResult, newTam); // p1 + p3 + p6
        subtract(bResult, p2, c22, newTam); // c22 = p1 + p3 - p2 + p6
 
        // Grouping the results obtained in a single matrix:
        for (i = 0; i < newTam ; i++) {
            for (j = 0 ; j < newTam ; j++) {
                c[i][j] = c11[i][j];
                c[i][j + newTam] = c12[i][j];
                c[i + newTam][j] = c21[i][j];
                c[i + newTam][j + newTam] = c22[i][j];
            }
        }
 
        // deallocating memory (free):
        a11 = free_real_matrix(a11, newTam);
        a12 = free_real_matrix(a12, newTam);
        a21 = free_real_matrix(a21, newTam);
        a22 = free_real_matrix(a22, newTam);
 
        b11 = free_real_matrix(b11, newTam);
        b12 = free_real_matrix(b12, newTam);
        b21 = free_real_matrix(b21, newTam);
        b22 = free_real_matrix(b22, newTam);
 
        c11 = free_real_matrix(c11, newTam);
        c12 = free_real_matrix(c12, newTam);
        c21 = free_real_matrix(c21, newTam);
        c22 = free_real_matrix(c22, newTam);
 
        p1 = free_real_matrix(p1, newTam);
        p2 = free_real_matrix(p2, newTam);
        p3 = free_real_matrix(p3, newTam);
        p4 = free_real_matrix(p4, newTam);
        p5 = free_real_matrix(p5, newTam);
        p6 = free_real_matrix(p6, newTam);
        p7 = free_real_matrix(p7, newTam);
        aResult = free_real_matrix(aResult, newTam);
        bResult = free_real_matrix(bResult, newTam);
 
} // end of Strassen function
 
/*------------------------------------------------------------------------------*/
// function to sum two matrices
void sum(double **a, double **b, double **result, int tam) {
 
    int i, j;
 
    for (i = 0; i < tam; i++) {
        for (j = 0; j < tam; j++) {
            result[i][j] = a[i][j] + b[i][j];
        }
    }
}
 
/*------------------------------------------------------------------------------*/
// function to subtract two matrices
void subtract(double **a, double **b, double **result, int tam) {
 
    int i, j;
 
    for (i = 0; i < tam; i++) {
        for (j = 0; j < tam; j++) {
            result[i][j] = a[i][j] - b[i][j];
        }
    }   
}
 
/*------------------------------------------------------------------------------*/
// This function allocates the matrix using malloc, and initializes it. If the variable random is passed
// as zero, it initializes the matrix with zero, if it's passed as 1, it initializes the matrix with random
// values. If it is passed with any other int value (like -1 for example) the matrix is initialized with no
// values in it. The variable tam defines the length of the matrix.
double **allocate_real_matrix(int tam, int random) {
 
    int i, j, n = tam, m = tam;
    double **v, a;         // pointer to the vector
 
    // allocates one vector of vectors (matrix)
    v = (double**) malloc(n * sizeof(double*));
 
    if (v == NULL) {
        printf ("** Error in matrix allocation: insufficient memory **");
        return (NULL);
    }
 
    // allocates each row of the matrix
    for (i = 0; i < n; i++) {
        v[i] = (double*) malloc(m * sizeof(double));
 
        if (v[i] == NULL) {
            printf ("** Error: Insufficient memory **");
            free_real_matrix(v, n);
            return (NULL);
        }
 
        // initializes the matrix with zeros
        if (random == 0) {
            for (j = 0; j < m; j++)
                v[i][j] = 0.0;
        }
 
        // initializes the matrix with random values between 0 and 10
        else {
            if (random == 1) {
                for (j = 0; j < m; j++) {
                    a = rand();
                    v[i][j] = (a - (int)a) * 10;
                }
            }
        }
    }
 
    return (v);     // returns the pointer to the vector. 
}
 
/*------------------------------------------------------------------------------*/
// This function unallocated the matrix (frees memory)
double **free_real_matrix(double **v, int tam) {
 
    int i;
 
    if (v == NULL) {
        return (NULL);
    }
 
    for (i = 0; i < tam; i++) { 
        if (v[i]) { 
            free(v[i]); // frees a row of the matrix
            v[i] = NULL;
        } 
    } 
 
    free(v);         // frees the pointer /
    v = NULL;
 
    return (NULL);   //returns a null pointer /
}