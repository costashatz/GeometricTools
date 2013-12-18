#include <iostream>
#include "MatrixAlgorithms.h"
using namespace std;

//SIMPLE EXAMPLE

int main(int argc, char *argv[])
{
	Matrix<3,3> A = Matrix<3,3>(3.,-6.,-3,2.,0.,6.,-4.,7.,4.);
	Vector<3> B = Vector<3>(-3.,-22.,3.);
	Vector<3> x = solveGauss(A,B);
	cout<<"Solving System: \n";
	A.Print();
	cout<<" = \n";
	B.Print();
	cout<<"\nSOLUTION: ";
	x.Print();
	cout<<"------------------------------------\n\n";
	Matrix<4,4> a = Matrix<4,4>(2,1,1,0,4,3,3,1,8,7,9,5,6,7,9,8);
	Vector<4> b = Vector<4>(1.,2.,3,4);
	cout<<"SOLUTION with LU: ";
	Vector<4> y = solveLU(a,b);
	y.Print();
	Vector<4> y2 = solveGauss(a,b);
	cout<<"SOLUTION with Gauss: ";
	y2.Print();
	cout<<"------------------------------------\n\n";
	return 0;
}